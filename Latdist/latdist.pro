;+
; NAME:
;
;   call_latdist
;
; PURPOSE:
;   Loops over csv indices and calls latdist which extrapolates the fit data to the centrifugal equator
;
; CALLING SEQUENCE:
;
;   call_latdist
;
; INPUTS: For this procedure, only keyword arguments are used
; 
; KEYWORS:
;   plot_extrapolation- plots input densities vs returned densities
;   save_data - saves the original densities as well as extrapolated data(including loaction) to a csv file
;   latdist - Walks up and down the field line model(which is hardcoded in).
;
;
; OUTPUTS:
;   CSV file containing spacecraft location in sIII right handed followed by the densities determined at that location. After which
;   the location of the centrifugal equator is given followed by the densities extrapolated to that point. 
;
; PROCEDURE:
;   Unpacks information about where the centrifugal equator is as well as the information about spacecraft position, densities and temperatures.
;   Once the data is unpacked, latdistVoy is called which is responsible for the actual calculation of the density values at the centrifugal equator.
;
; EXAMPLE:
; 
; 
; Created by Kaleb Bodisch

pro call_latdist, plot_extrapolation = plot_extrapolation, save_data = save_data, Field_Map = Field_Map

  degrad=!PI/180.

  ; get values of spacecraft location as well as cent equator location
  Field_line_model = read_csv('New_ephemeris_Field_Line_Map.csv', header = header)
  ;this flag tells you whether you are at 1-the start of a field line 2-the centrifugal equaotor of a field line 3-the spacecraft location along a given field line
  ;or 4-a location on a field line where the centrifugual equator location and spacecraft location are the same
  flag = Field_line_model.field1
  r_lines = Field_line_model.field2
  theta_lines = Field_line_model.field3
  phi_lines = Field_line_model.field4
  field_start_points = where(flag eq 1)
  
  number_lines = n_elements(field_start_points)
  
  theta_lines[where(theta_lines lt 90.)] = 90 - theta_lines[where(theta_lines lt 90.)]
  theta_lines[where(theta_lines ge 90.)] = -(theta_lines[where(theta_lines ge 90.)] -90)

  ;spacecraft measurement locations
  r_sc = r_lines[where(flag eq 3 or flag eq 4)]
  theta_sc = theta_lines[where(flag eq 3 or flag eq 4)]
  phi_sc = phi_lines[where(flag eq 3 or flag eq 4)]

  ;location of centrifugal equator using VIP4+CAN
  r_ce = r_lines[where(flag eq 2 or flag eq 4)]
  theta_ce = theta_lines[where(flag eq 2 or flag eq 4)]
  phi_ce = phi_lines[where(flag eq 2 or flag eq 4)]


  ;space_craft position for input into latdist
  scpos = transpose([[r_SC], [theta_sc]])

  ;co-latitude to latitude for input into latdist
  CEpos = transpose([[r_ce], [theta_ce]])

  ;order in which densities will be read into Latdist
  A = [16,16,32,32,32,1,23,16]
  Z = [1,2,1,2,3,1,1,1]
  
  
 
;this was for creating centrifugal equator profiles
;  CSV_Data = READ_CSV("V1_Data_CSV.csv", HEADER = CSVHeader)
;  ;get velocity in phi and organize densities/temperatures in the approriate order
;  v_phi = csv_data.field10
;  nion = [[csv_data.field15], [csv_data.field17], [csv_data.field19], [csv_data.field21], [csv_data.field23], [csv_data.field13], [csv_data.field25], [csv_data.field29]]   ;densities of species
;  tpar = [[csv_data.field16], [csv_data.field18], [csv_data.field20], [csv_data.field22], [csv_data.field24], [csv_data.field14], [csv_data.field26], [csv_data.field30]]

  
;make matrices to hold ion densities and temperatures
nion = make_array(8, number_lines)
tpar = make_array(8, number_lines)
  
  rho_ce = r_ce*cos(theta_ce*!dPI/180d)
  rho_sc = r_sc*cos(theta_sc*!dPI/180d)
  ;array for vphi
  v_phi = fltarr(number_lines)
  
  ;vphi profile produced from viper analysis
  x1 = where(rho_ce lt 8.53)
  x2 = where((rho_ce ge 8.53) and (rho_ce lt 17.005))
  x3 = where(rho_ce ge 17.005)
  v_phi[x1] =  rho_ce[x1]*12.572        ;assume corotation
  v_phi[x2] = 1.907*(rho_ce[x2]^2) + (-37.79*rho_ce[x2])+290.83286
  v_phi[x3] = .9666*rho_ce[x3] + 183.24
  
  
  fcorot = v_phi/(12.572*transpose(rho_ce))

  anis = 0d       ;isotropic
  BeB  = 0d       
  
  Aion = make_array(8, number_lines)
  Zion = make_array(8, number_lines)
  
  for l = 0, 7 do begin
    Aion[l, *] = A[l]
    Zion[l, *] = Z[l]  
  endfor
  
  ; calculate density and temperature from curves produced from viper reanlysis
  for k = 0, number_lines -1 do begin
    hot_flag = 0
    for i = 0, 7 do begin
      if ((Aion[i] eq 16) and (Zion[i] eq 2)) then begin
        a = 76.3 & b = -6.73& c = .086 & tpar[i, k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if ((Aion[i] eq 32) and (Zion[i] eq 1)) then begin
        a = 163  & b = -6.81& c = .169 & tpar[i, k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if ((Aion[i] eq 32) and (Zion[i] eq 2)) then begin
        a = 538. & b = -6.74& c = .598 & tpar[i, k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if ((Aion[i] eq 32) and (Zion[i] eq 3)) then begin
        a = 90.7 & b = -6.21& c = .165 & tpar[i, k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if ((Aion[i] eq 1) and (Zion[i] eq 1)) then begin
        a =50.6   & b = -5.31& c = .212 & tpar[i, k] = 94.1*(rho_ce[k]/6)^(.14)
      endif
      if ((Aion[i] eq 23) and (Zion[i] eq 1)) then begin
        a = 97.2 & b = -6.75& c = .106 & tpar[i, k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if ((Aion[i] eq 16) and (Zion[i] eq 1) and (hot_flag eq 1)) then begin
        a = 134.& b = -4.63& c = 1.057  & tpar[i,k] = 362*(rho_ce[k]/6)^(.91)
      endif
      if ((Aion[i] eq 16) and (Zion[i] eq 1) and (hot_flag eq 0)) then begin
        a = 592.& b = -7.36& c = .368 &hot_flag =1 & tpar[i,k] = 79.3*(rho_ce[k]/6)^(.714)
      endif
      if (rho_ce[k] lt 15.2) then nion[i,k] = a*(rho_ce[k]/6.)^b else nion[i,k] = c*(1987.*((rho_ce[k]/6.)^(-8.2))  + 14.*((rho_ce[k]/6.)^(-3.2)) + .05*((rho_ce[k]/6.)^(-.65)))
      if ((rho_ce[k] gt 16.25) and (Aion[i] eq 1) and (Zion[i] eq 1)) then nion[i,k] = rho_ce[k]*(-23./2375.) + (387./950.) ;linear curve for protons where fit curve doesnt work. 
    endfor
  endfor
  
  stop
  
  ;make sure these are cast as double
  Aion   = Double(Aion)
  Zion   = Double(Zion)
;  nion = transpose(nion)
;  tpar = transpose(tpar)
  v_phi = v_phi[0:(n_elements(r_sc)-1)]
  nion = nion[*,0:(n_elements(r_sc)-1)]
  tpar = tpar[*,0:(n_elements(r_sc)-1)]
  
  ;for centrifugal calculations
  output_struct = make_Array(9, n_elements(r_sc))
  single_field_point_results = make_array(9)
  total_point_number =long(0) 
   RETURN_ARR = FLTARR(11,n_elements(flag))       ;create a matrix to store values for r, z, each species density and total_charge_density
  for line = 0, number_lines-1 do begin  ;iterate over number of field lines
  ;for an individual field line, loop over number of points on that line
    if line lt number_lines-1 then field_points = [[r_lines[field_start_points[line]:(field_start_points[line+1])-1]],[theta_lines[field_start_points[line]:(field_start_points[line+1])-1]]] $
    else field_points = [[r_lines[field_start_points[line]:-1]],[theta_lines[field_start_points[line]:-1]]]
    field_points = transpose(field_points)
    
    print, line
    ns     = N_ELEMENTS(no)            ;number of species
    if keyword_set(field_map) then n_loops = n_elements(field_points[0,*]) else n_loops =1
    ;nloops will be 1 if we are just doing centrifugal equator else it its number of field line points
    for loop = 0, n_loops-1 do begin
      ;remove nans from nion and tpar
      tpar[where(finite(nion[*,line],/nan) eq 1,/null), loop] = 0
      nion[where(finite(nion[*,line],/nan) eq 1,/null), loop] = 0 


      ;only pass latdist values gt 0 else you will get nans in calculations
      ind = where(nion[*,line] gt 0, /null)
      tpar_lat = tpar[ind,line]
      Aion_lat = Aion[ind,line]
      Zion_lat = Zion[ind,line]
      nion_lat = nion[ind,line]
      ;include total electron density at given location in the output
      ind = [ind,8]

      ;if we are walking along field lines then save output into giant matrix else one array for each ceq point
      if keyword_set(field_map) then nion_CE = latdistVoy(cepos[*,line], field_points[*,loop], nion_lat, tpar_lat, Aion_lat, Zion_lat, fcorot[line], anis, BeB) else begin
        nion_CE = latdistVoy(scpos[*,line], CEpos[*,line], nion_lat, tpar_lat, Aion_lat, Zion_lat, fcorot[line], anis, BeB)
        output_struct[ind, line]= nion_ce
      endelse
      

      if keyword_set(Field_Map) then begin
        single_field_point_results[ind] = nion_ce
        RETURN_ARR(0,total_point_number)    = field_points(0,loop)
        RETURN_ARR(1,total_point_number)    = field_points(0,loop)*sin(field_points(1,loop)*!pi/180.) ; store z value
        RETURN_ARR(2:-1,total_point_number) = single_field_point_results
        total_point_number += 1
      endif
    endfor

  endfor  
  
  A = [16,16,32,32,32,1,23,16]
  Z = [1,2,1,2,3,1,1,1]

  if keyword_set(field_map) then write_csv, 'field_distribution_ephemeris_final.csv', return_arr
  
  if keyword_set(save_data) then  write_csv, 'CE_Data.csv', transpose([[r_sc],[theta_sc],[phi_sc],[transpose(nion[0,*])],[transpose(nion[1,*])],[transpose(nion[2,*])],[transpose(nion[3,*])],[transpose(nion[4,*])],[transpose(nion[5,*])],[transpose(nion[6,*])],[transpose(nion[7,*])], [r_ce], [theta_ce], [phi_sc], [transpose(output_struct[0,*])], [transpose(output_struct[1,*])], [transpose(output_struct[2,*])], [transpose(output_struct[3,*])], [transpose(output_struct[4,*])], [transpose(output_struct[5,*])], [transpose(output_struct[6,*])], [transpose(output_struct[7,*])], [transpose(output_struct[8,*])]]), Header=["r_sc","theta_sc","phi_sc","Oplus","Odouble","Splus","Sdouble","Stripple","Hplus","Naplus","Oplus_hot","r_ce","theta_ce","phi_ce","Oplus_ce","Odouble_ce","Splus_ce","Sdouble_ce","Stripple_ce","Hplus_ce","Naplus_ce","Oplus_hot_ce", "ntotal_ce"]

  if keyword_set(plot_extrapolation) then begin
    species = ["oplus", "oudble","splus","sdouble","stripple","Hplus","Na+", "oplus_hot"]
    for i =0, 7 do begin
        p = plot(r_sc[0:599], transpose(nion[i,[0:599]]),title = species[i],symbol = '.',linestyle ='',color = 'black', xrange=[4,40],xtitle = 'radius(x,y,z)',ytitle ="n/cc", yrange =[.01, 10000],/ylog ,/xlog)
        p = plot(r_sc[0:599], transpose(output_struct[i,[0:599]]),title = species[i],symbol = '.',linestyle ='',color = 'red', /overplot)
       endfor
  endif
end


FUNCTION latdistVoy, scpos, pos, no, tpar, mass, char, fcorot, anis, BeB, LATDIST=LATDIST, ANISOTROPIC=ANISOTROPIC, Field_Map = Field_Map

  ; Inputs
  ;  scpos - R(jupiter radii) and theta values of the spacecraft in jovigraphic SIII
  ;  pos - R(jupiter radii) and theta values of the Centrifugal equator in jovigraphic SIII
  ;  no - densities ion species in /cc
  ;  tpar - temperature of ions(in same order as no)
  ;  mass - A of given species(in same order as no)
  ;  char - charge state of each ion(in same order as no)
  ;  fcorot- fraction of corotation of plasma
  ;  anis - 0 for isotropic, how anisotropic plasma is
  ;   
  ; Keyword Arguments
  ;
  ; Outputs
  ;   Returns an array of extrapolated densities for individual species as well as total charge density
  ; 
  ;Created by Fran Bagenal
  ;Modified by Kaleb Bodisch
  ; 
  
  degrad=!PI/180.


;    ;new calculations from fran
    IF scpos(0) LT 12. THEN te   = 4.6*(scpos(0)/6)^3.4 ELSE te   = 50d
    IF scpos(0) LT 11.795734 THEN teh  = 35*(scpos(0)/6)^4.2 ELSE teh  = 600d
    IF scpos(0) LT 09.696 THEN nhnc = .001*(scpos(0)/6)^8 ELSE nhnc = 0.1d
   

  ; Check the lower limit boundary condition on the first Oxygen species - No limit on the Teh because nh -> 0 /cc
  IF tpar(0) LE 3 THEN BEGIN
    te   = tpar(0)
    nhnc = 0d
  ENDIF

  ;variable for number of species
  ns     = N_ELEMENTS(no)
  ;the following section will override individual densities for the densities determined by power laws
  ;if keyword_set(Field_Map) then begin
;   rho = scpos[0]*cos(scpos[1]*!pi/180.)
;   hot_flag = 0
;
;   for i = 0, ns - 1 do begin
;    if (mass[i] eq 16 and char[i] eq 2) then a = 76.3 & b = -6.73& c = .086 & tpar[i] = 79.3*(rho/6)^(.714)
;    if (mass[i] eq 32 and char[i] eq 1) then a = 163  & b = -6.81& c = .169 & tpar[i] = 79.3*(rho/6)^(.714)
;    if (mass[i] eq 32 and char[i] eq 2) then a = 538. & b = -6.74& c = .598 & tpar[i] = 79.3*(rho/6)^(.714)
;    if (mass[i] eq 32 and char[i] eq 3) then a = 90.7 & b = -6.21& c = .165 & tpar[i] = 79.3*(rho/6)^(.714)
;    if (mass[i] eq 1 and char[i] eq 1) then a =50.6   & b = -5.31& c = .212 & tpar[i] = 94.1*(rho/6)^(.14)
;    if (mass[i] eq 23 and char[i] eq 1 ) then a = 97.2& b = -6.75& c = .106 & tpar[i] = 79.3*(rho/6)^(.714)
;    if (mass[i] eq 16 and char[i] eq 1 and hot_flag eq 1) then a = 134.& b = -4.63& c = 1.057  & tpar[i] = 362*(rho/6)^(.91)
;    if (mass[i] eq 16 and char[i] eq 1 and hot_flag eq 0) then a = 592.& b = -7.36& c = .368 &hot_flag =1 & tpar[i] = 79.3*(rho/6)^(.714)
;    if rho[i] lt 15.2 then no[i] = a*(rho/6.)^b else no[i] = c*(1987.*((rho/6.)^(-8.2))  + 14.*((rho/6.)^(-3.2)) + .05*((rho/6.)^(-.65)))
;    if (rho gt 20. and mass[i] eq 1 and char[i] eq 1) then no[i] = .03
;    stop
;   endfor
; stop



  ; Create arrays of correct size
  ni     = FLTARR(ns+1)
  ni(ns) = TOTAL(no[*]*char[*]) ; Total electron density due to quasineutrality
  FOR i = 0, ns-1 DO ni(i) = no(i) ; Code feeds in density of individual species already - no need to calculate their density from fractions
 

    ngrid = nmax_Voy(ns,pos,scpos,ni,tpar,te,mass,char,fcorot,teh,nhnc)
    ; ngrid = nanismax_Voy(ns,pos,scpos,ni,tpar,te,mass,char,fcorot,anis,BeB,teh,nhnc)
  RETURN, ngrid
END

FUNCTION nmax_Voy, ns, pos, scpos, ni, tpar, te, mass, char, fcorot, teh, nhnc

  ; MAXWELLIAN DISTRIBUTIONS - ISOTROPIC
  ; includes effect of sub-corotational motion in reducing the
  ; centrifugal force by fcorot^2
  ; INPUT:
  ; ns = No. of ion species
  ; pos=position - where you want the densities calculated
  ; pos(0)=r
  ; pos(1)=lat
  ; scpos=position - where reference densities are measured/provided
  ; scpos(0)=r
  ; scpos(1)=lat
  ; ni(ns)=densities measured/provided at s/c (or reference location)
  ;        where ni(ns+1)=electron density
  ; tpar(ns)=parallel temperatures of ions
  ; te=electron temperature
  ; mass(ns)mass of each ion species
  ; char(ns)charge of each ion species
  ; Vphi/Vcorotation
  ; OUTPUT:
  ;  n(ns+1)= output densities of ns ion species plus n(ns+1)=electron density


  degrad=!PI/180.
  P=0.0
  n=fltarr(ns+1)

  ;CENTRIFUGAL FORCE ; Term 2 in final equation
  c2=pos(0)*pos(0)*cos(pos(1)*degrad)*cos(pos(1)*degrad)
  c1=scpos(0)*scpos(0)*cos(scpos(1)*degrad)*cos(scpos(1)*degrad)
  Cent=c1-c2
  ;P = -Cent

  ; START ITERATING
  it=20                   ;typical tests converge in less than 20 iterations.

  for k=0,it do begin
    f=0.
    df=0.

    ; ION SPECIES


    for i=0,ns-1 do begin
      e1=char(i)*P/tpar(i)
      ; jupiter e2=fcorot*fcorot*0.825*mass(i)*Cent/tpar(i)
      e2=fcorot*fcorot*0.825*mass(i)*Cent/tpar(i)
      ee=-e1-e2
      n(i)=ni(i)*exp(ee)
      f=f+char(i)*n(i)
      df=df-char(i)*char(i)*n(i)/tpar(i)
    endfor

    ; ELECTRONS
    e1=-1.*P/te
    e1h=e1*te/teh

    ; jupiter e2=0.825/1836.*Cent/Te
    e2=fcorot*fcorot*0.825/1836.*Cent/Te
    e2h=e2*te/teh

    ee=-e1-e2
    eeh=-e1h-e2h
    nc=ni(ns)*exp(ee)*(1.-nhnc)
    nh=ni(ns)*exp(eeh)*nhnc
    n(ns)=nc+nh

    ; Find Net Charge of ions + electrons
    f=f-n(ns)
    df=df-n(ns)/te

    P=P-f/df

    if abs(F) lt 0.001 then begin
      return,n
    endif

  endfor
  return,n
end
