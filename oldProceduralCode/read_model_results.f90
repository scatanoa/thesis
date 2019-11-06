subroutine read_float_arr(fileunit,path, record, N_x, vect, Res) 
	!Variables de entrada
	integer, intent(in) :: fileunit, record, N_x
	character*500, intent(in) :: path
	!Variables de salida
	real, intent(out) :: vect(N_x)
	integer, intent(out) :: Res
!	print *, 'initially, Res= ', Res   
	!Lectura 
	open(fileunit,file=path,form='unformatted',status='old',access='direct',&
		& RECL=4*N_x)
		read(fileunit,rec=record,iostat=Res) vect
		if (Res.ne.0) print *, 'Error: Se ha tratado de leer un valor fuera del rango'
	close(fileunit)
end subroutine

