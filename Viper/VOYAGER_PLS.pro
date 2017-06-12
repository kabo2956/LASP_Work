FUNCTION VOYAGER_PLS, Structure, POD, cupss, MAKE_PLOT=MAKE_PLOT, UNCERTAINTIES=UNCERTAINTIES, Structure_Unc=Structure_Unc, $
                      N_CHI_PTS=N_CHI_PTS, N_CHI_SIG=N_CHI_SIG, COMPUTE_CHI_SPACE=COMPUTE_CHI_SPACE
  ; Inputs
  ; Structure - A structure containing information on how to handle all parameters in the fit
  ; POD - Directory name for where to save the plots - Creates an output folder if necessary
  ; PLOT - Keyword argument to produce a plot of the data
  ; UNCERTAINTIES - Keyword argument to calculate uncertainties for the given parameters
  ; N_CHI_PTS         - Keyword argument for number of points to use (Creates a square array of that size)
  ; N_CHI_SIG         - Keyword argument for number of sigma to deviate from best fit
  ; Structure_Unc - Keyword argument with uncertainties of structures


  ; Output
  ; Return_Value - Returns a structure that is of the same size as source expression
  ;                Values will have changed to either better fit parameters or uncertainties in the parameters 
  
  ; Unpack Structure for relevant analysis
  ;
  ;Written by: Kaleb Bodisch & Logan Dougherty
  
  
  Year  = Structure.Year
  DOY   = Structure.DOY
  Hour  = Structure.Hour
  Minute = Structure.Minute
  Second = Structure.Second
  Spacecraft = Structure.Spacecraft
  IF Structure.L_or_M_Mode EQ 1 THEN n_channels = 128 ELSE n_channels = 16
  
  ; Define Forward Functions so IDL doesn't crash
  FORWARD_FUNCTION PLS_Data_Read
  FORWARD_FUNCTION Generate_PLS_Current
  
  ; Read in the Data for the timestamp and spacecraft given
  Dum_Struct = PLS_Data_Read(Structure) ; Works similar to Current_Plot_3 I believe (LOGAN 1/5/15)
  ; Combine data into one array, organized as a 4 x Number of Channels array
  PLS_Data_Transp = [Dum_Struct.current_array_A, Dum_Struct.current_array_B, Dum_Struct.current_array_C, Dum_Struct.current_array_D]
  ; Combine data into one array of 1 x (4 x Number of Channels) - It goes A, B, C, and then D cup
  PLS_Data = [TRANSPOSE(PLS_Data_Transp[0,*]),TRANSPOSE(PLS_Data_Transp[1,*]), $
              TRANSPOSE(PLS_Data_Transp[2,*]),TRANSPOSE(PLS_Data_Transp[3,*])]
  
  
  
  ; Calculate the errors in the Data
  ; Taken from Voyager Memo LOGAN 1/5/15 - PUT ACTUAL MEMO IN HERE THAT WE HAVE
  PLS_Error = SQRT(1d5 + 0.25d-4*PLS_Data*PLS_Data)
  
  ; If fitting procedure is turned on then fit the data
  IF Structure.Fit EQ 1 THEN BEGIN
    ; Remove data points that are not being fit (sets them to NaN, which MPFIT ignores)
    ; Create dummy array called Current_Holder
    Current_Holder = MAKE_ARRAY(n_channels*4,/DOUBLE)
    Current_Holder[*] = !VALUES.D_NAN        ; Initially populate with NaN values and then return true values

    ; Use Good_Channels function to find gchn 
    ; LOGAN 1/5/15 - MAKE SURE TO CHANGE MODE TO ACCOUNT FOR 0 AND 1 FOR L AND M
    cup = 0
    gchn = Good_Channels(long(Structure.Channels_CupA), cup, Structure.L_or_M_Mode)
    IF structure.Channels_CupA GT 0 THEN Current_Holder[gchn[0]-1 + 0*n_channels:gchn[1]-1 + 0*n_channels] $
                                                        = PLS_Data_Transp[0, gchn[0]-1:gchn[1]-1] 
    cup = 1
    gchn = Good_Channels(long(Structure.Channels_CupB), cup, Structure.L_or_M_Mode)
    IF structure.Channels_CupB GT 0 THEN Current_Holder[gchn[0]-1 + 1*n_channels:gchn[1]-1 + 1*n_channels] $
                                                        = PLS_Data_Transp[1, gchn[0]-1:gchn[1]-1]
    cup = 2
    gchn = Good_Channels(long(Structure.Channels_CupC), cup, Structure.L_or_M_Mode)
    IF structure.Channels_CupC GT 0 THEN Current_Holder[gchn[0]-1 + 2*n_channels:gchn[1]-1 + 2*n_channels] $
                                                        = PLS_Data_Transp[2, gchn[0]-1:gchn[1]-1]
    cup = 3
    gchn = Good_Channels(long(Structure.Channels_CupD), cup, Structure.L_or_M_Mode)
    IF structure.Channels_CupD GT 0 THEN Current_Holder[gchn[0]-1 + 3*n_channels:gchn[1]-1 + 3*n_channels] $
                                                        = PLS_Data_Transp[3, gchn[0]-1:gchn[1]-1]
    
    ; Set low values to NaNs so that they are not fit in the process
    ; Currently using a cutoff of 100 fA in the data to remove points where the instrument
    ; did not correctly measure the data - This is a user determined cutoff at the moment and there
    ; may be a better choice - Logan 3/1/16
    
    ind_remove = WHERE(Current_Holder LE 100, nind_remove)
    IF nind_remove GE 1 THEN Current_Holder(ind_remove) = !VALUES.D_NAN
    
    ; These points will also be removed from the calculation for the Chi-Space and the uncertainties
    
    ; Save to array for MPFIT to use
    PLS_Data_Fit = Current_Holder
    
    ; Create New Structure to house where parameters are located
    P_Locations = Structure
    FOR j = 0, N_TAGS(Structure) - 1 DO P_Locations.(j) = 0d ; Initially set all to zero
    
    ; Determine what velocities are active
    IF Structure.Vary_V1 EQ 1d THEN P_Locations.V1 = 1d
    IF Structure.Vary_V2 EQ 1d THEN P_Locations.V2 = 1d
    IF Structure.Vary_V3 EQ 1d THEN P_Locations.V3 = 1d
    
    ; Determine if both density and flowspeeds are allowed to vary
    IF Structure.Vary_Density EQ 1d THEN Vary_Density = 1d ELSE Vary_Density = 0d
    IF Structure.Vary_Temperature EQ 1d THEN Vary_Temperature = 1d
    IF Structure.Common_Temperature GT 0d AND Structure.Vary_Temperature NE 0 THEN BEGIN 
      Common_Temp = Structure.Common_Temperature
      P_Locations.Common_Temperature = 1d
    ENDIF
    
    ; Determine how many species there are
    StringCompare = STRCMP('Species_', TAG_NAMES(Structure), 3, /FOLD_CASE)
    ind = WHERE(StringCompare EQ 1, nind)
    n_species = nind/5
    Arguments = MAKE_ARRAY(n_species) ; Stores the arguments
    
    n_tied = 0 ; Counter for number of tied density parameters
    FOR j = 1, n_species DO BEGIN
      str1 = STRCOMPRESS(STRING('Species_', j), /REMOVE_ALL)
      StringCompare = STRCMP(str1, TAG_NAMES(Structure), 10, /FOLD_CASE) ; Extract type of parameter
      ; If it is species 10 or greater then there is one more letter in the variable name for IDL to check
      IF n_species GE 10 THEN StringCompare = STRCMP(str1, TAG_NAMES(Structure), 11, /FOLD_CASE)
      ind = WHERE(StringCompare EQ 1, nind)
      IF nind LT 1 THEN MESSAGE, 'Error: No species were provided for analysis. Quitting program.'
      Arguments(j-1) = Structure.(ind)

      
      ; If 0, then it is not in parameter space, either density or temperature (will follow its own temperature if common
      ; temperature isn't provided)
      IF Structure.(ind) EQ 0 THEN CONTINUE ELSE $
        
      ; If 1, then it is not in parameter space, and explicitly follows its own temperature
      IF Structure.(ind) EQ 1 THEN CONTINUE ELSE $
        
      ; If 2, then it is in parameter space, it will follow common temperature, or its own if common isn't provided
      IF Structure.(ind) EQ 2 THEN BEGIN
        ; Put density into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
        
        ; Check for common temperature tag - If not using common temperature, then include temperature in parameter space
        IF Structure.Common_Temperature EQ 0d THEN BEGIN
          str3 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
          str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind3 = WHERE(str3cmp EQ 1)
          P_Locations.(ind3) = 1d
        ENDIF 
      ENDIF ELSE IF Structure.(ind) EQ 3 THEN BEGIN
        ; Put density into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d

        ; Put temperature into parameter space
        str3 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
        str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind3 = WHERE(str3cmp EQ 1)
        P_Locations.(ind3) = 1d
      ENDIF ELSE IF (Structure.(ind) EQ 4) AND (Vary_Density EQ 1) THEN BEGIN
        ; Only the first tied density is a parameter in parameter space
        IF n_tied GT 0d THEN CONTINUE
        n_tied = n_tied + 1
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
      ENDIF ELSE IF Structure.(ind) EQ 5 THEN BEGIN
        ; Put temperature into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
        
        ; Put density into the parameter space
        IF n_tied GT 0d THEN CONTINUE
        IF Vary_Density EQ 1 THEN BEGIN
          n_tied = n_tied + 1
          str3 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
          str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind3 = WHERE(str3cmp EQ 1)
          P_Locations.(ind3) = 1d
        ENDIF
      ENDIF
    ENDFOR
    
    
    ; Save information to common blocks so it can be accessed inside of function that MPFIT calls
    Save_Structures, Struct_In_A=Structure, Struct_In_B=P_Locations
    
    ; Determine initial parameters
    dum_arr = MAKE_ARRAY(N_TAGS(Structure))
    FOR k = 0, N_TAGS(Structure) - 1 DO dum_arr(k) = P_locations.(k)
    ind = WHERE(dum_arr EQ 1)
    P = MAKE_ARRAY(N_ELEMENTS(ind), /DOUBLE)
    
    ; Generate parinfo structure for MPFIT to limit densities and temperatures to be above 0
    parinfo = REPLICATE({limited:[0,0],limits:[0.D,0]}, N_ELEMENTS(ind))

    FOR j = 0, N_ELEMENTS(ind) - 1 DO BEGIN
      P[j] = DOUBLE(Structure.(ind[j]))
      ; Check where density and temperature tags are set so that MPFIT can set limits and prevent values
      ; from going negative

      str1_check = '*_N*'
      str2_check = '*_T*'
      
      str1_check_match = STRMATCH((TAG_NAMES(Structure))[ind[j]], str1_check, /FOLD_CASE)
      str2_check_match = STRMATCH((TAG_NAMES(Structure))[ind[j]], str2_check, /FOLD_CASE)

      ; Check that it is a temperature or it is a density - Constrain it to above 0
      IF str1_check_match EQ 1 OR str2_check_match EQ 1 THEN BEGIN
        parinfo[j].limited[0] = 1
        parinfo[j].limits[0] = FLOAT(0.01)
      ENDIF
      
      ; Constrain temperatures to be under a certain value
      ;IF str2_check_match EQ 1 THEN BEGIN
      ;  parinfo[j].limited[1] = 1
      ;  parinfo[j].limits[1]  = FLOAT(300) ; 300 eV limit on temperatures in the fitting procedure
      ;ENDIF
    ENDFOR

    ; Need to set the parinfo command so that it doesn't set densities and temperatures to zero
    IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
      
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), ':', ULONG(Structure.DOY), ':', ULONG(Structure.Hour), ':', $
                                    ULONG(Structure.Minute), '::', ULONG(Structure.Second), $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02)'), /REMOVE_ALL)
      PRINT, STRCOMPRESS(STRING('Calculating parameter uncertainties for ', str_date))
      ; Determine what data points should be ignored in Chi-Square calculation
      index = WHERE(FINITE(PLS_Data_Fit) EQ 0, nindex) ; Indexes of calculation to avoid
      IF nindex EQ 0 THEN index = -9d ; Set as a filler number so that it will account for all data if none is bad
      Delta = 0.01d ; A 1 percent variation if the relative keyword is set
      Model_Name = 'Generate_PLS_Current'
      
      PLS_Unc = PLS_UNCERTAINTIES(Structure, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, Delta, /RELATIVE)
      Structure_Unc = Structure ; Set a similar structure for the uncertainties
      FOR j = 0, N_ELEMENTS(ind) - 1 DO Structure_Unc.(ind[j]) = PLS_Unc[j] ; Save uncertainties to structure
      
    ENDIF ELSE IF KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), ':', ULONG(Structure.DOY), ':', ULONG(Structure.Hour), ':', $
                                    ULONG(Structure.Minute), '::', ULONG(Structure.Second), $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02)'), /REMOVE_ALL)
      PRINT, STRCOMPRESS(STRING('Calculating Chi-Square space for ', str_date))
      ; Use underscores for text filename
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), '_', ULONG(Structure.DOY), '_', ULONG(Structure.Hour), '_', $
                                    ULONG(Structure.Minute), '_', ULONG(Structure.Second), '_CHI_ARRAYS.txt', $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)
      
      ; Run through the Chi-Square calculations
      N_Chi_Steps = N_Chi_Pts
      N_Sigma     = N_Chi_Sig
      Model_Name = 'Generate_PLS_Current'
      index = WHERE(FINITE(PLS_Data_Fit) EQ 0, nindex) ; Indexes of calculation to avoid
      IF nindex EQ 0 THEN index = -9d ; Set as a filler number so that it will account for all data if none is bad

      Structure_Hold = Structure ; Hold a copy of the structure before computing the chi arrays, just in case IDL overwrites it
      FORWARD_FUNCTION CHI_ARRAYS
      Chis = CHI_ARRAYS(Structure, Structure_Unc, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, N_Chi_Steps, N_Sigma, str_date)    
      Structure = Structure_Hold ; Return normal structure after computing chi arrays
    ENDIF ELSE BEGIN
      ; Run through MPFIT and calculate the best parameters

      ; Fit with MPFIT - Use all float for this part of the calculation - LOGAN 1/6/15 MPFITFUN does not seem to work with DOUBLE
      ; precision. Parameters are not changed how they should be for some reason. The flow velocities seem like they can change, yet
      ; everything else is completely stuck and does not vary. Simply changing to FLOAT precision allows for MPFIT to fit the
      ; parameters as it should, so it is not a limit that is being put on, but rather a limitation of MPFIT. It may be due to the
      ; range of the data for the PLS instrument, but it shouldn't be a coding issue on the end of the VIPER routines

      Best_P = MPFITFUN('Fit_Parameters', FLOAT(MAKE_ARRAY(4*n_channels, /DOUBLE)), FLOAT(PLS_Data_Fit), $
                        FLOAT(PLS_Error), FLOAT(P), /NOCATCH, /NAN, PARINFO=parinfo, MAXITER=Structure.Iterations) ; QUIET Keyword to suppress MPFIT from printing
      
      ; Now put the best parameters back into the structure
      FOR j = 0, N_ELEMENTS(ind) - 1 DO Structure.(ind[j]) = Best_P[j]
      ; Exit fitting routine
    ENDELSE
  ENDIF
  
  IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
    ; Generate final result of model
    PLS_Currents = Generate_PLS_Current(Structure)
    
    ; Sum up currents to get a total current curve
    IF Structure.L_or_M_Mode EQ 1 THEN n_channels = 128 ELSE n_channels = 16 ; Number of channels per cup
    Total_Current = MAKE_ARRAY(n_channels*4, /DOUBLE)

    ; Sum up currents from different ion species
    Modified_Array = TOTAL(PLS_Currents, 2) ; Reduce it to a 4 x Number of Channel Array, summing together the number of species
    
    ; Set up plot arguments
    xtitle = 'Channel Number'
    COMMON GENERATE
    IF (LDIST EQ "FALSE")  THEN ytitle = "Current (fA)"
    IF (LDIST EQ "TRUE") THEN ytitle = "Reduced Distribution Function (FemtoMhos)"

    ; Color titles for the plots
    colors = ['deep pink']
    ; Use Generate PLS Currents to get the plasma properties
    Plasma_Properties = Generate_PLS_Current(Structure, /RETURN_PLASMA_PROPERTIES)
    Masses       = Plasma_Properties[0,*]
    Charges      = Plasma_Properties[1,*]
    Densities    = Plasma_Properties[2,*]
    Temperatures = Plasma_Properties[3,*] 
    
    nO = 0 ; Count the number of Oxygen species - 1 species uses O+ as a proxy for hot species
    FOR nnn = 0, N_ELEMENTS(Masses) - 1 DO BEGIN
      A0 = Masses(nnn)
      Z0 = Charges(nnn)
      ; LOGAN 5/5 - Copying case statements from other section of code
      CASE A0 OF
        1:  colors = [colors, 'Black']
        16: CASE Z0 OF
              1: BEGIN
                   nO = nO + 1
                   IF nO EQ 1 THEN T0 = Temperatures(nnn)
                   T1 = Temperatures(nnn)
                   IF nO EQ 1 THEN colors = [colors, 'Dark Green']
                   IF nO EQ 2 THEN colors = [colors, 'Yellow Green']
                   IF T1 LT T0 THEN PRINT, "'Hot O+' species is at a lower temperature than primary O+ species. Both species will appear as a dark green color in the output plot."
                 END
              2: colors = [colors, 'Deep Sky Blue']
            ENDCASE
        23: CASE Z0 OF
              1: colors = [colors, 'Dark Orange']
            ENDCASE
        32: CASE Z0 OF
              3: colors = [colors, 'Red']
              2: colors = [colors, 'Maroon']
              1: colors = [colors, 'Dark Violet']
            ENDCASE
        64: CASE Z0 OF
              1: colors = [colors, 'Hot Pink']
            ENDCASE
        ELSE: colors = [colors, 'Gray']
      ENDCASE
    ENDFOR
    
    
    ; Positions for each faraday cup
    POSITIONS = [[.10,.60,.45,.90],[.60,.60,.95,.90], [.10,.20,.45,.50], [.60,.20,.95,.50]]

    ; Automatically Determine a range of currents (logarithmic scaling) that matches well with the data
    ;YRANGE = PLOT_RANGES_LOG(PLS_Data_Transp)
    YRANGE = [1d2,1d7] ; - Manual override
    ;YRANGE = [1d1,1d5]
    
    ; Turn off IDL errors for logarithmic plotting of values that may be zero
    !Except = 0
    
    ; Create dual error array so that array elements don't go below 0 - IDL will not plot negative values on a logarithmic scale
    PLS_Error_Array = MAKE_ARRAY(2, N_ELEMENTS(PLS_Error), /DOUBLE)
    
    ind_lt0 = WHERE(PLS_DATA - PLS_Error LE 0, nind_lt0)
    Dum_Arr = PLS_Error
    Dum_Arr[ind_lt0] = PLS_DATA[ind_lt0] - 1d-10 ; Set error values just above zero so that they appear in plot
    PLS_Error_Array(0,*) = Dum_Arr[*]
    PLS_Error_Array(1,*) = PLS_Error[*]
    ; Set errors to be in the proper format of 2 x 4 x number of channels
    Errors = MAKE_ARRAY(2,4,n_channels)
    FOR j = 0, 3 DO BEGIN
      Errors[*,j,*] = PLS_Error_Array[*, j*n_channels:(j+1)*n_channels-1]
    ENDFOR
    Cups = ['A','B','C','D']
    
    
    
    IF cupss eq 0 then begin
    ; Plot the Data
    ; Begin the loop through each cup
    FOR iii = 0, 3 DO BEGIN
      ; Values greater than 100 are for saving plots without them popping up
      IF Structure.Save_Plot GE 100 THEN BEGIN
        PLOT_BUFFER=1
      ENDIF
      dummy = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
                   xrange = [1,n_channels], POSITION = positions(*,iii), BUFFER=PLOT_BUFFER, yrange = YRANGE, Color = 'Deep Pink') 
      dummy_dat = PLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], /OVERPLOT, POSITION = positions(*,iii), thick=2)
      eplot = ERRORPLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
          ERRORBAR_COLOR = 'Orange', $
          ERRORBAR_CAPSIZE = 0.0, $
          errorbar_thick = 1, $
          /OVERPLOT)
      ; Calculate the number of species and then overplot them individually
      n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
      FOR jjj = 0, n_species-1 DO BEGIN
        pplot= PLOT(FINDGEN(n_channels)+1, PLS_Currents(iii,jjj,*), title = Cups[iii]+' Cup', xtitle = xtitle, color=colors[jjj+1],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
                     ytitle = ytitle, yrange = yrange, /ylog, /overplot)
      ENDFOR
      dummy_tot = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
                       xrange = [1,n_channels], POSITION = positions(*,iii), BUFFER=PLOT_BUFFER, yrange = YRANGE, Color = 'Deep Pink', /OVERPLOT, thick=1.5)
    ENDFOR
    endif
    
    IF cupss NE 0 THEN BEGIN
      positionss = [.10,.20,.95,.90]
      IF cupss EQ 1 THEN Cup = 'A'
      IF cupss EQ 2 THEN Cup = 'B'
      IF cupss EQ 3 THEN Cup = 'C'
      IF cupss EQ 4 THEN Cup = 'D'
                        
      

      ; Plot the Data
      ; Begin the loop through each cup
      FOR iii = cupss-1, cupss-1 DO BEGIN
        ; Values greater than 100 are for saving plots without them popping up
        IF Structure.Save_Plot GE 100 THEN BEGIN
          PLOT_BUFFER=1
        ENDIF
        dummy = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720],  $
          xrange = [1,n_channels], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE, Color = 'Deep Pink')
        dummy_dat = PLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], thick=2,/overplot)
        eplot = ERRORPLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
          ERRORBAR_COLOR = 'Orange', $
          ERRORBAR_CAPSIZE = 0.0, $
          errorbar_thick = 1, $
          /OVERPLOT)
        ; Calculate the number of species and then overplot them individually
        n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
        FOR jjj = 0, n_species-1 DO BEGIN
          pplot= PLOT(FINDGEN(n_channels)+1, PLS_Currents(iii,jjj,*), title = Cups[iii]+' Cup', xtitle = xtitle, color=colors[jjj+1],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
            ytitle = ytitle, yrange = yrange, /ylog, /overplot)
        ENDFOR
        dummy_tot = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720],  $
          xrange = [1,n_channels], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE, Color = 'Deep Pink', /OVERPLOT, thick=1.5)
      ENDFOR
    ENDIF

    ; Print the relevant information on the graphic
    
    ; Determine Response, Radial Distance, and Time of plot
    IF Structure.Response EQ 0 THEN Response_Name = 'CUPINT/DCPINT' ELSE $
    IF Structure.Response EQ 1 THEN Response_Name = 'LABCUR/LDCUR'
    
    IF Structure.Spacecraft EQ 1 THEN BEGIN
      COMMON VGR1Data
      SSEDR_Data = VGR1SSEDR
    ENDIF ELSE IF Structure.Spacecraft EQ 2 THEN BEGIN
      COMMON VGR2Data
      SSEDR_Data = VGR2SSEDR
    ENDIF
    
    ; Time in SSEDR file
    SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d
    
    ; Time of Input
    Time_Requested = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
    ; Time of Data
;    IF Structure.L_or_M_Mode EQ 0 THEN Data_Time = Dum_Struct.Ltimes_Array ELSE $
;    IF Structure.L_or_M_Mode EQ 1 THEN Data_Time = Dum_Struct.Mtimes_Array
    
    Data_Time = Dum_Struct.times_Array
    
    ; Create DOY, Hour, Minute, and Second from Time_Data
    ; LOGAN 1/7/15 - Check to see if this works
    Data_DOY  = FIX(Data_Time)
    Data_Hour = FIX((Data_Time - Data_DOY)*24d)
    Data_Min  = FIX(((Data_Time - Data_DOY)*24d - Data_Hour)*60d)
    Data_Sec  = FIX((((Data_Time - Data_DOY)*24d - Data_Hour)*60d - Data_Min)*60d)
    
    ; Find SSEDR index of data to find S/C location
    SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
    ; Get X, Y, and Z locations of spacecraft and compute radial distance (Rj)
    xSC_S3 = SSEDR_Data[6,SSEDR_ind] ; System III cartesian coordinates
    ySC_S3 = SSEDR_Data[7,SSEDR_ind]
    zSC_S3 = SSEDR_Data[8,SSEDR_ind]
    Rj = SQRT(xSC_S3*xSC_S3+ySC_S3*ySC_S3+zSC_S3*zSC_S3)/71492d

    ; Print information to graphic
    t = TEXT(0.5,0.98, 'Mode:' + Response_Name + '           Year:' + STRING(ULONG(Structure.Year)) + $
                       '           Date (DOY, Hour, Min, Second):' + $
                       STRING(ULONG(Data_DOY)) + STRING(ULONG(Data_Hour)) + STRING(ULONG(Data_Min)) + STRING(ULONG(Data_Sec)) + $
                       '           RJ:    ' + STRING(Rj, format='(F0.3)'),  font_size = 10, color = 'black', $
                       ALIGNMENT = .5, vertical_ALIGNMENT = .5, /normal)
    ; Print remaining argument tags
    u = TEXT(.02, .09, 'A (amu), Z (q):', /normal)
    u = TEXT(.02, .06, '$n (cm^{-3}):$',  /normal)
    u = TEXT(.02, .03, 'T (eV):',  /normal)

    ; Flowspeed
    q=text(.02, .12, '$Cyl Vel(V_r, V_{\phi}, V_z):$')  
    dum_names = [Structure.V1, Structure.V2, Structure.V3] 
    sep=.05
    FOR k = 0, 2 DO q=text((1+k)*sep*2d + 0.15,0.12, string(dum_names[k], format='(F0.2)')) 

    ; Use Generate PLS Currents to get the plasma properties  
    Plasma_Properties = Generate_PLS_Current(Structure, /RETURN_PLASMA_PROPERTIES)
    Masses       = Plasma_Properties[0,*]
    Charges      = Plasma_Properties[1,*]
    Densities    = Plasma_Properties[2,*]
    Temperatures = Plasma_Properties[3,*]
    
    ; Determine seperation based on the number of species
    seps = [0.29,0.26,0.23,0.20,0.17,0.15,0.12,0.10,0.08]
    seperation = seps[N_ELEMENTS(Masses)-1]
    
    FOR k = 0, N_ELEMENTS(Masses)-1 DO BEGIN
      v = TEXT((1+k)*seperation+.047,0.09, STRCOMPRESS(STRING(UINT(Masses(k)))+','+STRING(UINT(Charges(k)))), color=colors[k+1], /norm)
      v = TEXT((1+k)*seperation+.05,0.06, STRING(Densities(k),format='(F0.2)'), /norm)
      v = TEXT((1+k)*seperation+.05,0.03, STRING(Temperatures[k], format='(F0.2)'), /norm)
    ENDFOR
    
    ; Save a plot if requested
    IF Structure.Save_Plot GE 1 THEN BEGIN
      
      ; Create name for plot based on current UTC time
      YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
      DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
      IF ULONG(DayOfRun) LT 10 THEN BEGIN
        DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
      ENDIF
      MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
      HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
      MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
      SecondOfRun = STRMID(SYSTIME(/UTC),17,2)
      
      ; Name the file
      IF ((Structure.Save_Plot EQ 1) OR (Structure.Save_Plot EQ 101)) THEN $
      OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
                              HourOfRun, MinuteOfRun, '_', SecondOfRun, '.jpg', $
                              FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
      ELSE IF ((Structure.Save_Plot EQ 2) OR (Structure.Save_Plot EQ 102)) THEN $
      OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
                              HourOfRun, MinuteOfRun, '_', SecondOfRun, '.png', $
                              FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
      ELSE IF ((Structure.Save_Plot EQ 3) OR (Structure.Save_Plot EQ 103)) THEN $
      OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
                              HourOfRun, MinuteOfRun, '_', SecondOfRun, '.ps', $
                              FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
      ELSE PRINT, 'Unable to save selected image file type.'

      
      ; Save the file if of the 3 types - Can add different image types if desired at a later date
      IF ISA(OutputFilename) NE 0 THEN BEGIN
        ; Change directories for output of plots
        CD, CURRENT = RESTORE_DIR
        CD, POD
        dummy.SAVE, OutputFilename, Border = 10, RESOLUTION = 150
        ; Change back directories
        CD, RESTORE_DIR
      ENDIF
    ENDIF
    
    ; Plot the data, the errorbars, and the model
    ; If Structure.Save_Plot is on then just go ahead and save the plot automatically, with a predetermined name
    
    
  ENDIF
  
  ; If uncertainties are requested, return the uncertainty structure. Otherwise return the best parameters
  IF KEYWORD_SET(UNCERTAINTIES) THEN RETURN, Structure_Unc ELSE $
                                     RETURN, Structure

  
END

FUNCTION FIT_PARAMETERS, X, P, err=err
  ; This function is essential for MPFIT operations. MPFIT calls on it and the parameters are then
  ; passed into the structure accordingly. Once they are in the structure, the function calls on the
  ; model function to return the currents. This function sums up the currents from each species
  ; so that it can return just a total current to MPFIT and MPFIT can iteratively change the parameters
  ; 
  ; Inputs
  ; X - Array of X values to be used by MPFIT
  ; P - Array of parameter values to use
  ; err - Error values on the Data
  ; 
  ; Additional Inputs:
  ; These are passed in via a common block so that the function can have the proper
  ; structure format
  ; 
  ; Structure    - A structure with the values to use - has to be augmented with P that is passed in by MPFIT
  ; P_Locations  - A structure with the locations of the parameters so that they can easily be put into the structure
  ; 
  ; Outputs:
  ; Total_Current - A 512 (M-Mode) or 64 (L-Mode) array that has the total current
  
  ; Read in Structure and P_Locations from a common block that was saved before this function was called
  COMMON Structure_Block
  Structure   = varA
  P_Locations = varB
  ; Put parameters correctly into the Structure variable by using P_Locations
  dum_arr = MAKE_ARRAY(N_TAGS(Structure))
  FOR k = 0, N_TAGS(Structure) - 1 DO dum_arr(k) = P_locations.(k)
  ind_fit = WHERE(dum_arr EQ 1)
  FOR j = 0, N_ELEMENTS(ind_fit) - 1 DO Structure.(ind_fit[j]) = P[j]
  ; Use the PLS model to generate the current
  
  FORWARD_FUNCTION Generate_PLS_Current
  Total_Current = Generate_PLS_Current(Structure, /RETURN_TOTAL_CURRENT, /SET_FLOAT)
  
  ; Return Total_Current to MPFIT so it can find best fit parameters
  RETURN, Total_Current
END

FUNCTION Generate_PLS_Current, Structure, RETURN_TOTAL_CURRENT=RETURN_TOTAL_CURRENT, SET_FLOAT=SET_FLOAT, $
  RETURN_PLASMA_PROPERTIES=RETURN_PLASMA_PROPERTIES, calculate_chi = calculate_chi
  ; This function calculates the current given the structure input. It will read elements from the structure
  ; and use the given composition to produce the simulated current that the program will read. It depends on
  ; both the planet and the type of response (Maxwellian, Kappa, etc.) that is specified
  ;
  ; Currently works for the following planets:
  ; Jupiter
  ;
  ; Currently works with the following responses:
  ; CUPINT/DCPINT - Convected Maxwellian
  ; LABCUR/LDCUR  - Convected Maxwellian (accounts better for a warm plasma distribution)
  ;
  ; Inputs:
  ; Structure - A structure that holds all of the pertinent information for the model. Information about
  ;             contents of the structure can be found close to the top of VIPER.pro


  ; Determine Spacecraft and pull appropriate SSEDR information for trajectory and location
  ; All lookup tables are included in CONRD in common blocks
  IF Structure.Spacecraft EQ 1 THEN BEGIN
    COMMON VGR1DATA
  ENDIF ELSE BEGIN
    COMMON VGR2DATA
  ENDELSE


  ; Jupiter based analysis
  IF Structure.Planet_Number EQ 5 THEN BEGIN

    ; LOGAN 1/5/15
    ; Currently, both CUPINT and LABCUR use essentially the same inputs so they will both be
    ; set up the same. There may need to be an IF loop here to setup a Kappa distribution differently
    ; though likely not
    ; LOGAN 1/5/15

    ; Setup array for the number of species that current is being calculated for
    ; If the density of a species is 0, then it will not be accounted for because that just
    ; adds zero to the current
    ;
    ; Store arrays of the important arguments (A, Z, n, T, and how the ion is being processed)

    ; Calculate the number of ions
    StringCompare = STRCMP('Species_', TAG_NAMES(Structure), 3, /FOLD_CASE)
    ind = WHERE(StringCompare EQ 1, nind)
    n_species = nind/5

    ; Establish arrays to put elements in
    Arguments    = MAKE_ARRAY(n_species)
    Masses       = Arguments
    Charges      = Arguments
    Densities    = Arguments
    Temperatures = Arguments
    
    
    ; Fill all arrays
    FOR k = 0, n_species - 1 DO BEGIN
      str0 = STRCOMPRESS(STRING('Species_', k+1), /REMOVE_ALL)
      IF k+1 LT 10 THEN str0cmp = STRCMP(str0, TAG_NAMES(Structure), 10, /FOLD_CASE) ELSE $
                        str0cmp = STRCMP(str0, TAG_NAMES(Structure), 11, /FOLD_CASE)
      ind0 = WHERE(str0cmp EQ 1)
      Arguments(k) = Structure.(ind0)

      str1 = STRCOMPRESS(STRING('Species_', k+1, '_A'), /REMOVE_ALL)
      IF k+1 LT 10 THEN str1cmp = STRCMP(str1, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str1cmp = STRCMP(str1, TAG_NAMES(Structure), 13, /FOLD_CASE)
      ind1 = WHERE(str1cmp EQ 1)
      Masses(k) = Structure.(ind1)

      str2 = STRCOMPRESS(STRING('Species_', k+1, '_Z'), /REMOVE_ALL)
      IF k+1 LT 10 THEN str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 13, /FOLD_CASE)
      ind2 = WHERE(str2cmp EQ 1)
      Charges(k) = Structure.(ind2)

      str3 = STRCOMPRESS(STRING('Species_', k+1, '_N'), /REMOVE_ALL)
      IF k+1 LT 10 THEN str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str3cmp = STRCMP(str3, TAG_NAMES(Structure), 13, /FOLD_CASE)
      ind3 = WHERE(str3cmp EQ 1)
      Densities(k) = Structure.(ind3)

      str4 = STRCOMPRESS(STRING('Species_', k+1, '_T'), /REMOVE_ALL)
      IF k+1 LT 10 THEN str4cmp = STRCMP(str4, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str4cmp = STRCMP(str4, TAG_NAMES(Structure), 13, /FOLD_CASE)
      ind4 = WHERE(str4cmp EQ 1)
      Temperatures(k) = Structure.(ind4)
    ENDFOR

    ; Check densities - If a species has 0 density, then it will not have a contribution to the current
    ind_good = WHERE(Densities GT 0, nind_good)

    ; Augment arrays to only contain relevant ion species
    Arguments    = Arguments(ind_good)
    Masses       = Masses(ind_good)
    Charges      = Charges(ind_good)
    Densities    = Densities(ind_good)
    Temperatures = Temperatures(ind_good)

    ; Fill in temperature and density arrays
    n_tied = 0
    FOR k = 0, nind_good - 1 DO BEGIN
      ; Change temperatures to a common temperature if the parameters are indicated to have those values
      IF Structure.Common_Temperature GT 0d THEN BEGIN
        IF ((Arguments(k) EQ 0) OR (Arguments(k) EQ 2) OR (Arguments(k) EQ 4)) THEN $
          Temperatures(k) = Structure.Common_Temperature
      ENDIF

      ; Check if density values should be tied together and set them accordingly
      IF ((Arguments(k) EQ 4) OR (Arguments(k) EQ 5)) THEN BEGIN
        n_tied = n_tied + 1
        IF n_tied EQ 1 THEN BEGIN
          Master_Density = Densities(k)
        ENDIF ELSE BEGIN
          Densities(k) = Master_Density*Densities(k)
        ENDELSE
      ENDIF
    ENDFOR

    ; Do a quick check that no temperatures are set to zero. If they are, that ion species will be ignored as
    ; it will lead to non finite results
    ind_good2 = WHERE(Temperatures GT 0d, nind_good2)
    IF nind_good2 NE nind_good THEN BEGIN
      Arguments    = Arguments(ind_good2)
      Masses       = Masses(ind_good2)
      Charges      = Charges(ind_good2)
      Densities    = Densities(ind_good2)
      Temperatures = Temperatures(ind_good2)
    ENDIF

    ; Return plasma properties if requested - Useful for plotting procedure
    IF KEYWORD_SET(RETURN_PLASMA_PROPERTIES) THEN BEGIN
      new_arr = MAKE_ARRAY(4, nind_good2, /DOUBLE)
      new_arr[0,*] = Masses[*]
      new_arr[1,*] = Charges[*]
      new_arr[2,*] = Densities[*]
      new_arr[3,*] = Temperatures[*]
      RETURN, new_arr
    ENDIF


    IF Structure.L_or_M_Mode EQ 1 THEN n_channels = 128 ELSE n_channels = 16 ; Number of channels per cup
    ; Set up array to store current information
    ; First element is the Cup (A, B, C, D)
    ; Second element is the current/ion species
    ; Third element is the channel number of the PLS instrument
    PLS_Currents = MAKE_ARRAY(4, nind_good2, n_channels, /DOUBLE)


    ; Get some of the constants for the arguments of the CUPINT/DCPINT and LABCUR/LDCUR programs
    ; Contains the following variables:
    ;   SS    - Step size for numerical integration
    ;   LDIST - Logical flag to return current values ('FALSE') or reduced current distribution ('TRUE')
    COMMON Generate

    ; Now that all species should be correctly accounted for, begin calculating the currents
    FOR j = 0, 3 DO BEGIN ; Begin the loop over the Faraday cup
      Cups = ['A','B','C','D']

      ; Put flow speeds into a vector
      VelocityCylindrical = [Structure.V1, Structure.V2, Structure.V3]

      ; Calculate decimal date of the time being used
      Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d

      ; Calculate flowspeed into Faraday cup
      U = VGR_CupVelocity(VelocityCylindrical, Cups[j], Structure) ; FILLER - NEED THIS FUNCTION


      ; Convert the temperatures into a thermal speed (km/s)
      ;FORWARD_FUNCTION eV_ThermalSpeed ;
      ; LOGAN 1/6/15 - eV_ThermalSpeed isn't being recognized for some reason
      Thermals = eV_ThermalSpeed(Temperatures, Masses) ; FILLER - NEED THIS FUNCTION

      ; Begin loop over relevant species
      FOR k = 0, nind_good2 - 1 DO BEGIN

        ; Calculate response of the instrument
        ; CUPINT/DCPINT response
        IF Structure.Response EQ 0 THEN BEGIN
          
          ; Calculate currents for A, B, and C cups
          IF j NE 3 THEN BEGIN
            PLS_Currents(j,k,*) = CUPINT(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode,[1,n_channels],SS,LDIST)
            ; Calculate response for D cup
          ENDIF ELSE BEGIN
            PLS_Currents(j,k,*) = DCPINT(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode,[1,n_channels],SS,LDIST)
          ENDELSE

        ; LABCUR/LDCUR response
        ENDIF ELSE IF Structure.Response EQ 1 THEN BEGIN
          
          ; Calculate Currents for A, B, and C cups
          IF j NE 3 THEN BEGIN
            PLS_Currents(j,k,*) = LABCUR(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode, $
              [1,n_channels],LDIST,Cups[j])
            ; Calculate response for D cup
          ENDIF ELSE BEGIN
            PLS_Currents(j,k,*) = LDCUR(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode, $
              [1,n_channels],LDIST,Cups[j])
          ENDELSE
        
        ; CUPINT/DCPINT with a Kappa distribution addition
        ENDIF ELSE IF Structure.Response EQ 2 THEN BEGIN
          
          ; Not programmed to work
          PRINT, 'No Kappa distribution has been programmed yet. Returning currents of 0'
          PLS_Currents(j,k,*) = 0
          
        ENDIF ELSE BEGIN
          PRINT, 'Error: There is not a response function for this. Returning a current of zero for all channels'
          PLS_Currents(j,k,*) = 0
        ENDELSE

      ENDFOR
    ENDFOR


    ; Have the option to return the current already summed
    IF KEYWORD_SET(RETURN_TOTAL_CURRENT) THEN BEGIN
      ; It always goes in the order of Cup A, Cup B, Cup C, and then Cup D
      ; PLS_Currents is 4 x Number of Species x Number of Channels
      IF Structure.L_or_M_Mode EQ 1 THEN n_channels = 128 ELSE n_channels = 16 ; Number of channels per cup
      Total_Current = MAKE_ARRAY(n_channels*4, /DOUBLE)

      ; LOGAN 1/5/15 - CHECK THAT THIS SUM WORKS
      Modified_Array = TOTAL(PLS_Currents, 2) ; Reduce it to a 4 x Number of Channel Array, summing together the number of species
      ; Modify Total_Current to have either 512 (M-Mode) or 64 (L-Mode) successive channels
      Total_Current(0*n_channels:1*n_channels-1) = Modified_Array(0, 0:n_channels-1)
      Total_Current(1*n_channels:2*n_channels-1) = Modified_Array(1, 0:n_channels-1)
      Total_Current(2*n_channels:3*n_channels-1) = Modified_Array(2, 0:n_channels-1)
      Total_Current(3*n_channels:4*n_channels-1) = Modified_Array(3, 0:n_channels-1)
      IF KEYWORD_SET(SET_FLOAT) THEN RETURN, FLOAT(TOTAL_CURRENT) ELSE RETURN, Total_Current
    ENDIF ELSE BEGIN
      ; Return array of currents per individual species
      IF KEYWORD_SET(SET_FLOAT) THEN RETURN, FLOAT(PLS_Currents) ELSE RETURN, PLS_Currents
    ENDELSE

  ENDIF

END

FUNCTION PLS_DATA_READ, Structure
  ; For Jupiter 
  IF Structure.Planet_Number EQ 5 THEN BEGIN
    
    ; Create string for directory of text files, as well as the current directory
    CD, CURRENT = RESTORE_DIR
    STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Jupiter/PLS_Data'))
    
    ; Determine DOY and Time
    DOY = Structure.DOY
    Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
    
    ; Name of text file that contains pertinent data
    Data_Filename = STRCOMPRESS(STRING('day', LONG(DOY), '.txt'), /REMOVE_ALL)
    FORWARD_FUNCTION PLS_Current_Data
    
    
    RETURN, PLS_Current_Data(Structure, Data_Filename, Time, 'D', STR_DIR)
  ENDIF
END

FUNCTION PLS_Current_Data, Structure, Data_Filename, Time, Cup, STR_DIR
  FORWARD_FUNCTION read_file
  data = read_file(Data_Filename, CHANGE_DIRECTORY=STR_DIR)    ; read in the data from the file "filename"
  
  Mode = Structure.L_or_M_Mode
  ; make energy arrays for M modes, L modes, and the E modes:

  ; get array of energies using the conversionM function
  FORWARD_FUNCTION conversionM
  ; There are 128 channels for M modes - each corresponds to an energy bin
  MchannelsLow = findgen(128)             ; lower bound on bins (channel number 0-127)
  MchannelsHigh = 1. + findgen(128)       ; upper bound on bins (channel number 1-128)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  Menergies = ( conversionM(MchannelsLow) + conversionM(MchannelsHigh) ) / 2.   ; convert channel numbers to energy values for M modes

  ; get array of energies using the conversionLE2 function
  FORWARD_FUNCTION conversionLE2
  ; There are 16 channels for L modes
  LchannelsLow = findgen(16)              ; lower bound on bins (channel number 0-15)
  LchannelsHigh = 1. + findgen(16)        ; upper bound on bins (channel number 1-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  Lenergies = ( conversionLE2(LchannelsLow) + conversionLE2(LchannelsHigh) ) / 2.    ; convert channel numbers to energy values for L modes

  ; make energy arrays for both E modes using conversion functions
  FORWARD_FUNCTION conversionLE2
  FORWARD_FUNCTION conversionE1
  ; There are 16 channels for E1 modes
  E1channelsLow = findgen(16)             ; lower bound on bins (0-15)
  E1channelsHigh = 1. + findgen(16)       ; upper bound on bins (1-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  E1energies = ( conversionE1(E1channelsLow) + conversionE1(E1channelsHigh) ) / 2.  ; convert channel numbers to energy values for E1 modes
  ; There are 16 channels for E2 modes, but the first 4 are garbage, so use only the last 12
  E2channelsLow = 4. + findgen(12)        ; lower bound on bins (4-15)
  E2channelsHigh = 5. + findgen(12)       ; upper bound on bins (5-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  E2energies = ( conversionLE2(E2channelsLow) + conversionLE2(E2channelsHigh) ) / 2. ; convert channel numbers to energy values for E2 modes
  ; now, combine E1 and E2 energies into one big energy array
  Eenergies = fltarr(28)                   ; make one array for all energies
  Eenergies[0:15] = E1energies             ; first half of energy array is the E1 energy levels
  Eenergies[16:27] = E2energies            ; other half of energy array is the E2 energies
  
  
  FORWARD_FUNCTION convert_times
  ; use mode chosen by user to define the important quantities for the plots (currents, energies, and times)
  IF (mode EQ 1) THEN BEGIN                  ; if the user chose M modes
    modetimes = data.times.M                     ; important time string array
    datatimes = convert_times(data.times.M)      ; time strings converted to decimal day-of-year (to compare with user's chosen time)
    energies = Menergies                         ; important energy array
    Acurrents = data.cupA.M.toarray()            ; cup A currents are M modes
    Bcurrents = data.cupB.M.toarray()            ; cup B currents are M modes
    Ccurrents = data.cupC.M.toarray()            ; cup C currents are M modes
    Dcurrents = data.cupD.M.toarray()            ; cup D currents are M modes

    current__holder =make_array(3612)
    current_total= [Acurrents,Bcurrents,Ccurrents,Dcurrents]
    current_total = round(current_total)
    current__holder = current_total[uniq(current_total, sort(current_total))]
    current__array_size = (size(current__holder))[1]
    xvals=findgen(current__array_size)
    ;pp = plot(xvals, current__holder, /ylog, ytitle = 'unique values', xtitle = 'elements')
    ;ppl = plot(xvals, current__holder, ytitle = 'unique values', xtitle = 'elements')
    ;gradient =alog(current__holder[*])/xvals[*]
    ;print, gradient
    ;gradplot = plot(xvals, gradient, ytitle = 'gradient', xtitle = 'elements')

  ENDIF ELSE IF (mode EQ 0) THEN BEGIN       ; if the user chose L modes
    modetimes = data.times.L                     ; important time string array
    datatimes = convert_times(data.times.L)      ; time strings converted to decimal day-of-year (to compare with user's chosen time)
    energies = Lenergies                         ; important energy array
    Acurrents = data.cupA.L.toarray()            ; cup A currents are L modes
    Bcurrents = data.cupB.L.toarray()            ; cup B currents are L modes
    Ccurrents = data.cupC.L.toarray()            ; cup C currents are L modes
    Dcurrents = data.cupD.L.toarray()            ; cup D currents are L modes

    Around = round(Acurrents)
    Bround = round(Bcurrents)
    Cround = round(Ccurrents)
    Dround = round(Dcurrents)

    Acurrentsunique = [uniq(Around, sort(Around))] ; sorting arrays and removing duplicate values
    Bcurrentsunique = [uniq(Bround, sort(Bround))] ; sorting arrays and removing duplicate values
    Ccurrentsunique = [uniq(Cround, sort(Cround))] ; sorting arrays and removing duplicate values
    Dcurrentsunique = [uniq(Dround, sort(Dround))] ; sorting arrays and removing duplicate values
    current__holder =make_array(3612)
    current_total= [Acurrentsunique,Bcurrentsunique,Ccurrentsunique, Dcurrentsunique]


    current_total = round(current_total)
    current__holder = current_total

    current__array_size = (size(current__holder))[1]
    xvals=findgen(current__array_size)



  ENDIF ELSE IF (mode EQ 'E') THEN BEGIN       ; if the user chose E modes
    modetimes = data.times.E1                    ; time string array (going with E1 modes for this since the 2 E-modes will be combined)
    datatimes = convert_times(data.times.E1)     ; time strings converted to decimal day-of-year (to compare with user's chosen time)
    energies = Eenergies                         ; important energy array
    ; combine the currents
;    E1currents = data.Emodes.E1.toarray()        ; array for E1 currents
;    E2currents = data.Emodes.E2.toarray()        ; array for E2 currents
;    ; data may be missing for one mode and not for another, so ensure that the E2 mode to be plotted with each E1 mode
;    ;   is from a time close enough to the E1 mode that it must be from the same "scan"
;    t = n_elements(data.times.E1)                ; define total number of times to be number of times of E1 mode data
;    Ecurrents = dblarr(t, 28)                    ; make a 2D double array big enough to hold all E data (16 E1 channels + 12 E2 channels)
;    Ecurrents[*,0:15] = E1currents               ; first 16 elements of Ecurrents array correspond to the E1 data
;    E1times = convert_times(data.times.E1)      ; convert E1 times to decimal DOY to compare how close E2 times are
;    E2times = convert_times(data.times.E2)      ; convert E2 times to decimal DOY to compare how close they are to E1
;    ; There are 48 seconds between the time that an E1 mode is taken and the time an E2 mode is taken.
;    timespan = 4.8d1 / 8.64d4         ; this is the fraction of a day taken up by 48 seconds (the time between when E1 and E2 were taken)
;    ; match up E2 currents to where their corresponding E1 currents are in the Ecurrents array
;    FOR i = 0, t-1 DO BEGIN           ; for each E1 mode time
;      FOR j = 0, n_elements(E2times)-1 DO BEGIN     ; for each E2 mode time
;        ; if an E2 mode time falls within 48 seconds of an E1 mode time, the E2 currents correspond to those E1 currents
;        IF ( (E2times[j] LE (E1times[i] + timespan)) AND (E2times[j] GT (E1times[i] - timespan)) ) THEN BEGIN  ; if E2 data is near E1 data      ; removed 1/7 no need for e mode info
;          Ecurrents[i,16:27] = E2currents[j,*]       ; add these E2 currents to the Ecurrents array in the same spot as the E1 currents
;          ; if there's an E1 mode with no matching E2 mode, the E2 part of the spectrum will be considered 0 and not plotted
;        ENDIF       ; an E2 mode is within 48 seconds of an E1 mode
;      ENDFOR      ; each element of E2times
;    ENDFOR      ; each element of E1times
    ; now all E mode quantities are defined
  ENDIF      ; important quantities have been defined for every mode the user could have selected

  ; Get closest data time stamp
  times= closest(datatimes, time)

  new_stuff={current_array_A:Acurrents[times,*], current_array_B:Bcurrents[times,*], current_array_C:Ccurrents[times,*], current_array_D:Dcurrents[times,*], times_array:datatimes[times,*]}
  ;create a structure of current values and the corresponding time value for this measurement to be returned from the function
  return, new_stuff
END

FUNCTION read_file, Data_Filename, CHANGE_DIRECTORY=CHANGE_DIRECTORY
  IF KEYWORD_SET(CHANGE_DIRECTORY) THEN CD, CHANGE_DIRECTORY, CURRENT = RESTORE_DIR
  openr, lun, Data_Filename, /get_lun       ; open file
  ; make lists for each type of data that will be read in from the file
  E1currentArrays = list()        ; list for all arrays of E1 currents
  E2currentArrays = list()        ; list for all arrays of E2 currents
  cupA_LcurrentArrays = list()    ; list for all arrays of L currents in cup A
  cupA_McurrentArrays = list()    ; list for all arrays of M currents in cup A
  cupB_LcurrentArrays = list()    ; list for all arrays of L currents in cup B
  cupB_McurrentArrays = list()    ; list for all arrays of M currents in cup B
  cupC_LcurrentArrays = list()    ; list for all arrays of L currents in cup C
  cupC_McurrentArrays = list()    ; list for all arrays of M currents in cup C
  cupD_LcurrentArrays = list()    ; list for all arrays of L currents in cup D
  cupD_McurrentArrays = list()    ; list for all arrays of M currents in cup D
  E1times = list()        ; list for all times at which E1 data were taken
  E2times = list()        ; list for all times at which E2 data were taken
  Ltimes = list()         ; list for all times at which L data were taken
  Mtimes = list()         ; list for all times at which M data were taken

  ; Now actually start reading the file.
  WHILE ~EOF(lun) DO BEGIN
    ; read in time & garbage lines
    firstline = ""               ; predefine to be a string
    readf, lun, firstline        ; read in the first line of data
    ; This first line contains the time and status (which tells us the mode).
    ; Now there are 3 garbage lines:
    garbageLine0 = ""            ; predefine lines of garbage to be strings
    garbageLine1 = ""
    garbageLine2 = ""
    readf, lun, garbageLine0     ; read in garbage lines (will be ignored)
    readf, lun, garbageLine1
    readf, lun, garbageLine2
    ; Alter the firstline string to split into time and status parts
    newFirstLine = strmid(firstline, 16)       ; get rid of beginning unnecessary characters
    splitFirstLine = strsplit(newFirstLine, 'status', /extract)  ; separate time & status
    time = splitFirstLine[0]               ; time is the first part of the split string
    status = strmid(splitFirstLine[1], 2)  ; status is second part of split string, minus extra spaces
    ; read status to determine next steps
    mode = strmid(status, 0, 1)  ; This tells us the first character of the status word, which gives the mode
    ; Here's how this "mode vs status word" thing works:
    ;   If the first character of the status word is:             then the mode is:
    ;                                                 (a space)                    L
    ;                                               1                              L
    ;                                               2                              L
    ;                                               3                              L
    ;                                               4                              E1
    ;                                               5                              E1
    ;                                               6                              E1
    ;                                               7                              E1
    ;                                               8                              E2
    ;                                               9                              E2
    ;                                               A                              E2
    ;                                               B                              E2
    ;                                               C                              M
    ;                                               D                              M
    ;                                               E                              M
    ;                                               F                              M
    ; Now we read the data into different lists based on what the mode is (using an IF test).
    ; The IF test below basically does the following for each section of data:
    ; puts time in whatever times list it belongs in
    ; makes an array to read currents into
    ; reads currents into said array
    ; tacks that array onto the end of whatever list it belongs in
    ; for L or M modes, repeats for every cup
    IF ((mode EQ '1') OR (mode EQ '2') OR (mode EQ '3') OR (mode EQ " ")) THEN BEGIN     ; L modes
      Ltimes.add, time                     ; add the time string to the L-modes times list
      ; L modes have 4 cups that measure 16 currents, with each array of currents separated by garbage in the file.
      ; predefine the 4 current arrays (of doubles) and 3 garbage lines:
      dummyCurrentA = dblarr(16)
      garbage0 = ""
      dummyCurrentB = dblarr(16)
      garbage1 = ""
      dummyCurrentC = dblarr(16)
      garbage2 = ""
      dummyCurrentD = dblarr(16)
      ; now read in the data and assign it to the appropriate lists
      readf, lun, dummyCurrentA                 ; read in cup A currents
      readf, lun, garbage0                      ; will be ignored
      readf, lun, dummyCurrentB                 ; read in cup B currents
      readf, lun, garbage1                      ; will be ignored
      readf, lun, dummyCurrentC                 ; read in cup C currents
      readf, lun, garbage2                      ; will be ignored
      readf, lun, dummyCurrentD                 ; read in cup D currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 15 DO BEGIN
        IF (dummyCurrentA[i] LE 1.) THEN dummyCurrentA[i] = 1.
        IF (dummyCurrentB[i] LE 1.) THEN dummyCurrentB[i] = 1.
        IF (dummyCurrentC[i] LE 1.) THEN dummyCurrentC[i] = 1.
        IF (dummyCurrentD[i] LE 1.) THEN dummyCurrentD[i] = 1.
      ENDFOR   ; every current in each of the four cups
      cupA_LcurrentArrays.add, dummyCurrentA    ; add this array to the end of the list of cup A's L-mode currents
      cupB_LcurrentArrays.add, dummyCurrentB    ; add this array to the end of the list of cup B's L-mode currents
      cupC_LcurrentArrays.add, dummyCurrentC    ; add this array to the end of the list of cup C's L-mode currents
      cupD_LcurrentArrays.add, dummyCurrentD    ; add this array to the end of the list of cup D's L-mode currents
    ENDIF ELSE IF ((mode EQ '4') OR (mode EQ '5') OR (mode EQ '6') OR (mode EQ '7')) THEN BEGIN   ; E1 modes
      E1times.add, time                    ; add the time string to the E1-modes times list
      ; E1 modes have 1 cup that measures 16 currents
      dummyCurrent = dblarr(16)            ; predefine array of 16 currents (doubles)
      readf, lun, dummyCurrent             ; read in currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 15 DO BEGIN
        IF (dummyCurrent[i] LE 1.) THEN dummyCurrent[i] = 1.
      ENDFOR   ; every current in the array
      E1currentArrays.add, dummyCurrent    ; add this array to the end of the list of the E1 currents
    ENDIF ELSE IF ((mode EQ '8') OR (mode EQ '9') OR (mode EQ 'A') OR (mode EQ 'B')) THEN BEGIN      ; E2 modes
      E2times.add, time                    ; add the time string to the E2-modes times list
      ; E2 modes have 1 cup that measures 16 currents
      dummyCurrent = dblarr(16)            ; predefine array of 16 currents (doubles)
      readf, lun, dummyCurrent             ; read in currents
      ; the first 4 channels are garbage, so get rid of them
      newCurrent = dummyCurrent[4:*]
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 11 DO BEGIN
        IF (newCurrent[i] LE 1.) THEN newCurrent[i] = 1.
      ENDFOR   ; every current in the array
      E2currentArrays.add, newCurrent    ; add this array to the end of the list of the E2 currents
    ENDIF ELSE IF ((mode EQ 'C') OR (mode EQ 'D') OR (mode EQ 'E') OR (mode EQ 'F')) THEN BEGIN   ; M modes
      Mtimes.add, time                     ; add the time string to the M-modes times list
      ; M modes have 4 cups that measure 128 currents, with each array of currents separated by garbage in the file.
      ; predefine the 4 current arrays (of doubles) and 3 garbage lines:
      dummyCurrentA = dblarr(128)
      garbage0 = ""
      dummyCurrentB = dblarr(128)
      garbage1 = ""
      dummyCurrentC = dblarr(128)
      garbage2 = ""
      dummyCurrentD = dblarr(128)
      ; now read in the data and assign it to the appropriate lists
      readf, lun, dummyCurrentA                 ; read in cup A currents
      readf, lun, garbage0                      ; will be ignored
      readf, lun, dummyCurrentB                 ; read in cup B currents
      readf, lun, garbage1                      ; will be ignored
      readf, lun, dummyCurrentC                 ; read in cup C currents
      readf, lun, garbage2                      ; will be ignored
      readf, lun, dummyCurrentD                 ; read in cup D currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 127 DO BEGIN
        IF (dummyCurrentA[i] LE 1.) THEN dummyCurrentA[i] = 1.
        IF (dummyCurrentB[i] LE 1.) THEN dummyCurrentB[i] = 1.
        IF (dummyCurrentC[i] LE 1.) THEN dummyCurrentC[i] = 1.
        IF (dummyCurrentD[i] LE 1.) THEN dummyCurrentD[i] = 1.
      ENDFOR   ; every current in each of the four cups
      cupA_McurrentArrays.add, dummyCurrentA    ; add this array to the end of the list of cup A's M-mode currents
      cupB_McurrentArrays.add, dummyCurrentB    ; add this array to the end of the list of cup B's M-mode currents
      cupC_McurrentArrays.add, dummyCurrentC    ; add this array to the end of the list of cup C's M-mode currents
      cupD_McurrentArrays.add, dummyCurrentD    ; add this array to the end of the list of cup D's M-mode currents
    ENDIF        ; every mode has been accounted for
  ENDWHILE       ; we haven't encountered the end of the file

  ; convert the lists of times into string arrays for use with convert_times function later
  LtimesArray = Ltimes.toarray(type=STRING)
  MtimesArray = Mtimes.toarray(type=STRING)
  E1timesArray = E1times.toarray(type=STRING)
  E2timesArray = E2times.toarray(type=STRING)
  ; Now, make anonymous structures to hold the data.
  cupA = {L: cupA_LcurrentArrays, M: cupA_McurrentArrays}    ; structure to hold L & M mode data for cup A
  cupB = {L: cupB_LcurrentArrays, M: cupB_McurrentArrays}    ; structure to hold L & M mode data for cup B
  cupC = {L: cupC_LcurrentArrays, M: cupC_McurrentArrays}    ; structure to hold L & M mode data for cup C
  cupD = {L: cupD_LcurrentArrays, M: cupD_McurrentArrays}    ; structure to hold L & M mode data for cup D
  Emodes = {E1: E1currentArrays, E2: E2currentArrays}        ; structure to hold E1 & E2 mode data
  times = {L: LtimesArray, M: MtimesArray, E1: E1timesArray, E2: E2timesArray}   ; structure to hold all time arrays
  ; now make one big anonymous structure to hold the other structures, so all info can be returned at once
  data = {cupA:cupA, cupB:cupB, cupC:cupC, cupD:cupD, Emodes:Emodes, times:times}
  free_lun, lun
  IF KEYWORD_SET(CHANGE_DIRECTORY) THEN CD, RESTORE_DIR
  return, data       ; return the structure full of data
END         ; read_file

FUNCTION read_info, filename

  openr, lun, filename, /get_lun           ; open the file

  ; define a set of lists into which to sort the information from the file
  times = list()                ; list for date/time information
  longitudes = list()           ; list for Jupiter longitudes along Voyager's trajectory
  latitudes = list()            ; list for Jupiter latitudes along Voyager's trajectory
  magneticLatitudes = list()    ; list for Jupiter's latitudes with respect to its magnetic poles along Voyager's trajectory
  localTimes = list()           ; list for local times on Jupiter along Voyager's trajectory
  radialDistances = list()      ; list for Voyager's radial distance from Jupiter (units of Jovian radii)
  IoPhases = list()             ; list for phases of Io at the times given

  ; now read in the data and sort it
  WHILE ~EOF(lun) DO BEGIN           ; until the end of the file is encountered
    line = ""              ; predefine line to be a string
    readf, lun, line       ; read one full line
    linePieces = strsplit( strcompress(line), " ", /extract)    ; split the line into substrings where there are spaces


    ; first five substrings are time: year, day-of-year, hour, minute, second
    ; combine these first five substrings into one "time" string
    time = linePieces[0] + " " + linePieces[1] + " " + linePieces[2] + " " + linePieces[3] + " " + linePieces[4]
    times.add, time              ; add this time string to the times list
    ; longitude is the next quantity in the array of substrings (linePieces)
    longitude = double( linePieces[5] )            ; convert the substring to a double and label it as longitude
    longitudes.add, longitude                      ; add this longitude to the longitudes list

    ; latitude is the next quantity in the array of substrings (linePieces)
    latitude = double( linePieces[6] )             ; convert the substring to a double and label it as latitude
    latitudes.add, latitude                        ; add this latitude to the latitudes list

    ; magnetic latitude is the next quantity in the array of substrings (linePieces)
    magneticLatitude = double( linePieces[7] )     ; convert the substring to a double and label it as magneticLatitude
    magneticLatitudes.add, magneticLatitude        ; add this magneticLatitude to the magneticLatitudes list

    ; local time is the next quantity in the array of substrings (linePieces)
    localTime = double( linePieces[8] )            ; convert the substring to a double and label it as localTime
    localTimes.add, localTime                      ; add this localTime to the localTimes list

    ; radial distance is the next quantity in the array of substrings (linePieces)
    radialDistance = double( linePieces[9] )       ; convert the substring to a double and label it as radialDistance
    radialDistances.add, radialDistance            ; add this radialDistance to the radialDistances list

    ; Io phase is the last quantity in the array of substrings (linePieces)
    IoPhase = double( linePieces[10] )             ; convert the substring to a double and label it as IoPhase
    IoPhases.add, IoPhase                          ; add this IoPhase to the IoPhases list

  ENDWHILE        ; there are still lines to be read in the file

  ; convert the lists to arrays, now that they're full of all the info in the file
  timesArray = times.ToArray(type='STRING', /no_copy)               ; convert list to a string array for the times
  longitudesArray = longitudes.ToArray(type='DOUBLE', /no_copy)     ; convert list to a double array for longitudes
  latitudesArray = latitudes.ToArray(type='DOUBLE', /no_copy)       ; convert list to a double array for latitudes
  magneticLatitudesArray = magneticLatitudes.ToArray(type='DOUBLE', /no_copy)  ; convert list to a double array for mag. latitudes
  localTimesArray = localTimes.ToArray(type='DOUBLE', /no_copy)     ; convert list to a double array for local times
  radialDistancesArray = radialDistances.ToArray(type='DOUBLE', /no_copy)  ; convert list to a double array for radial distances
  IoPhasesArray = IoPhases.ToArray(type='DOUBLE', /no_copy)         ; convert list to a double array for Io phases

  ; put each of these arrays into a structure, VoyagerInfo, as a field
  VoyagerInfo = {times:timesArray, longitudes:longitudesArray, latitudes:latitudesArray, magneticLatitudes:magneticLatitudesArray, $
    localTimes:localTimesArray, radialDistances:radialDistancesArray, IoPhases:IoPhasesArray}

  return, VoyagerInfo         ; return the structure filled with the info from the file (infoFileName)

END         ; read_info function

;------------------------------------------------------------------------
; conversion functions: conversionM, conversionLE2, conversionE1
;------------------------------------------------------------------------
;
;------------------------------------------------------------------------
;
; conversionM
;   function that converts an array of channels into an array of energies
;   (to be used with M-Mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 128 channels ranging from 1 to 128)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray
;            (steps in energyArray are logarithmic between 10eV and 5950eV)
;
;------------------------------------------------------------------------
FUNCTION conversionM, channelArray
  energyArray = 60. * 10.^(channelArray / 64.) - 50.    ; conversion formula for M modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionM

;
;------------------------------------------------------------------------
;
; conversionLE2
;   function that converts an array of channels into an array of energies
;   (to be used with L-mode or E2-mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 16 channels ranging from 1 to 16)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray
;            (steps in energyArray are logarithmic between 10eV and 5950eV)
;
;------------------------------------------------------------------------
FUNCTION conversionLE2, channelArray
  energyArray = 60. * 10.^(channelArray / 8.) - 50.    ; conversion formula for L and E2 modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionLE2

;
;------------------------------------------------------------------------
;
; conversionE1
;   function that converts an array of channels into an array of energies
;   (to be used with E1-mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 16 channels ranging from 1 to 16)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray140eV)
;
;------------------------------------------------------------------------
FUNCTION conversionE1, channelArray
  energyArray = 60. * 10.^(channelArray / 32.) - 50.    ; conversion formula for E1 modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionE1

;------------------------------------------------------------------------
; convert_times function
;
;  function to convert an array of strings representing times to an
;    array of doubles (decimal day of year)
;  to be used with data from Voyager
;
;  INPUT: timesStringArray - an array of strings that each give a time
;           format of string is: year, day of year, hour, minute, second
;           #### ### ## ## ##.####
;
;  OUTPUT: decimalDOYarray, the same array of times but converted
;            from strings into doubles with format of decimal DOY
;            (ie, day.fraction-of-day, ###.#######)
;
;  Note: because this is data from Voyager, the year is always 1979,
;           so the year will be ignored (not converted)
;
;------------------------------------------------------------------------
FUNCTION convert_times, timesStringArray
  decimalDOYarray = dblarr( n_elements(timesStringArray) )  ; make an array of doubles that's as big as the string array
  FOR time = 0, n_elements(timesStringArray)-1 DO BEGIN     ; for every element (every time) in the string array
    compressedTime = strcompress(timesStringArray[time])      ; get rid of extra spaces in the string
    timePieces = strsplit(compressedTime, " ", /extract)      ; split up the string into substrings wherever there are spaces
    ; Now timePieces is an array of substrings: [year, DOY, hour, minute, second]
    ; we know which element of timePieces corresponds to each kind of time unit, so convert each element into a double
    year = double(timePieces[0])               ; convert the first element of timePieces, year, into a double
    day = double(timePieces[1])                ; convert the DOY into a double
    hour = double(timePieces[2])               ; convert the hour into a double
    minute = double(timePieces[3])             ; convert the minute into a double
    second = double(timePieces[4])             ; convert the second into a double
    ; turn these into a value for decimal DOY by converting each part to fractions of a day, then adding them all together:
    ;    day / (1 day per day)
    ; +  hour / (24 hours per day)
    ; +  minute / (60 minutes per hour) / (24 hours per day)
    ; +  second / (60 seconds per minute) / (60 minutes per hour) / (24 hours per day)
    ;  = decimal DOY
    ; then fill in the corresponding element of decimalDOYarray with the result
    decimalDOYarray[time] = day + hour/24. + minute/(60.*24.) + second/(3600.*24.)
  ENDFOR         ; every time in timesStringArray, so we fill up decimalDOYarray
  return, decimalDOYarray     ; return the new array of times converted to doubles of decimal DOY
END     ; convert_times


