      program ASCComplex
c form perfect geometry of dimer if two monomers fall in distance cutoff
      implicit none

      integer num_SubUnit_template
      parameter (num_SubUnit_template=12)
      integer npart
      parameter (npart=8)
      integer mol_tot_num
      parameter (mol_tot_num=50)
      integer mol_type_num
      parameter (mol_type_num=1)
      integer interface_type_num
      parameter (interface_type_num=2)
      real*8 cell_range_x
      parameter (cell_range_x=100.0)
      real*8 cell_range_y
      parameter (cell_range_y=100.0)
      real*8 cell_range_z
      parameter (cell_range_z=100.0)
      integer simu_step
      parameter (simu_step=50000000)
      real*8 time_step
      parameter (time_step=0.1)
      real*8 distance_step
      parameter (distance_step=2.0)
      real*8 mol_D
      parameter (mol_D=1.0)
      real*8 mol_rot_D
      parameter (mol_rot_D=2.0)
      real*8 pai
      parameter (pai=3.1415926)  
      real*8 radius
      parameter (radius=1.25)
      real*8 tolerance_r
      parameter (tolerance_r=0.1)
      integer num_trajec
      parameter (num_trajec=20)
      real*8 bond_dist_cutoff
      parameter (bond_dist_cutoff=0.5)
      integer output_flag_rec,output_flag_str
      parameter (output_flag_rec=1,output_flag_str=1)   
      integer trj_output_freq
      parameter (trj_output_freq=10000)
      integer rec_output_freq
      parameter (rec_output_freq=1000)
      real*8 bond_thetapd_1,bond_thetapd_cutoff_1
      parameter (bond_thetapd_1=180.0,bond_thetapd_cutoff_1=30.0)
      real*8 bond_thetaot_1,bond_thetaot_cutoff_1
      parameter (bond_thetaot_1=60.0,bond_thetaot_cutoff_1=30.0)
      real*8 bond_thetapd_2,bond_thetapd_cutoff_2
      parameter (bond_thetapd_2=180.0,bond_thetapd_cutoff_2=30.0)
      real*8 bond_thetaot_2,bond_thetaot_cutoff_2
      parameter (bond_thetaot_2=0.0,bond_thetaot_cutoff_2=30.0)	  
      real*8 Flx_TlRg
      parameter (Flx_TlRg=0.25)! Unit: nm
      real*8 Flx_RtRg
      parameter (Flx_RtRg=10.0)! unit: degree	  
cc
      real*8 cm_SubUnit_x(num_SubUnit_template)
      real*8 cm_SubUnit_y(num_SubUnit_template)
      real*8 cm_SubUnit_z(num_SubUnit_template)
      integer i,j,index,i2,j2,k,k2,n_t
      real*8 theta,height
      real*8 dist
      real*8 res_x_c(npart),res_y_c(npart),res_z_c(npart)
      real*8 x_axis_bg,y_axis_bg,z_axis_bg
      real*8 x_axis_ed,y_axis_ed,z_axis_ed
      real*8 x_bf_rot,y_bf_rot,z_bf_rot
      real*8 x_af_rot,y_af_rot,z_af_rot
      real*8 res_x(npart,num_SubUnit_template)
      real*8 res_y(npart,num_SubUnit_template)
      real*8 res_z(npart,num_SubUnit_template)
      real*8 temp_i,temp_j,temp_k
      real*8 phi,psai,phai
      real*8 t(3,3)
      real*8 mol_x_0(npart,mol_tot_num)
      real*8 mol_y_0(npart,mol_tot_num)
      real*8 mol_z_0(npart,mol_tot_num)	  
      real*8 mol_x(npart,mol_tot_num)
      real*8 mol_y(npart,mol_tot_num)
      real*8 mol_z(npart,mol_tot_num)
      real*8 mol_x0(npart,mol_tot_num)
      real*8 mol_y0(npart,mol_tot_num)
      real*8 mol_z0(npart,mol_tot_num)
      real*8 mol_x0_new(npart,mol_tot_num)
      real*8 mol_y0_new(npart,mol_tot_num)
      real*8 mol_z0_new(npart,mol_tot_num)
      real*8 mol_x0_new0(npart,mol_tot_num)
      real*8 mol_y0_new0(npart,mol_tot_num)
      real*8 mol_z0_new0(npart,mol_tot_num)	  
      real*8 mol_x_new(npart,mol_tot_num)
      real*8 mol_y_new(npart,mol_tot_num)
      real*8 mol_z_new(npart,mol_tot_num)
      real*8 mol_x_new0(npart,mol_tot_num)
      real*8 mol_y_new0(npart,mol_tot_num)
      real*8 mol_z_new0(npart,mol_tot_num)
      integer mol_status(mol_tot_num)
      integer mol_status_new(mol_tot_num)
      integer site_status(npart,mol_tot_num)
      integer site_status_new(npart,mol_tot_num)
      integer mol_type_index(mol_tot_num)
      integer mol_index
      real*8 displace_amp
      real*8 initial_simu_time,current_simu_time
      integer mc_time_step
      integer iteration_mole_step
      real*8 PB_x,PB_y,PB_z
      integer collision_flag
      real*8 Prob_Diff,Prob_Ass,Prob_Diss
      integer complex_num,complex_num_new
      integer complex_mol_idx(mol_tot_num*npart,2)
      integer complex_site_idx(mol_tot_num*npart,2)
      integer complex_mol_idx_new(mol_tot_num*npart,2)
      integer complex_site_idx_new(mol_tot_num*npart,2)
      integer atompr(2,1000)
      real*8 xa(1000),ya(1000),za(1000),xb(1000),yb(1000),zb(1000)
      integer npair,npair_a,npair_b
      real*8 RMSD
      integer selected_complex
      integer selected_mol1
      integer selected_mol2
      integer selected_site1
      integer selected_site2
      integer align_num_i
      integer align_num_j
      integer align_ind
      integer target_ind_i2
      integer target_ind_j2
      real*8 Ass_Rate_max(mol_type_num,mol_type_num,interface_type_num)
      real*8 Diss_Rate_max(mol_type_num,mol_type_num,interface_type_num)
      real*8 Diss_Rate_Value
      real*8 Ass_Rate_Value
      integer oligomer_num
      integer oligomer_subunit_num(mol_tot_num)
      integer oligomer_subunit_index(mol_tot_num,mol_tot_num)
      integer oligomer_state(mol_tot_num)
      integer contact_matrix(mol_tot_num,mol_tot_num)
      integer color_flag(mol_tot_num)
      integer oligomer_map(mol_tot_num,mol_tot_num)
      integer count_flag
      real*8 point_x(3),point_y(3),point_z(3)
      real*8 theta_pd,theta_ot
      real*8 Flx_Rg_x,Flx_Rg_y,Flx_Rg_z
	  
      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5

      r3=5.0   
      r4=5.0
      r5=5.0    
	  
cccccccccccccccccccccccccccccccccccccccccccccccc
c    Read Kinetic Parameters
cccccccccccccccccccccccccccccccccccccccccccccccc

      open (unit=10,file=
     &     'CageKineticParameter_34.dat',
     &     status='old')
      do i=1,mol_type_num
         do j=1,mol_type_num
            do k=1,interface_type_num
               read(10,1100) Ass_Rate_max(i,j,k),Diss_Rate_max(i,j,k)
            enddo
         enddo
      enddo
      close(10)
 1100 format(9X,2F9.5)

ccccccccccccccccccccccccccccccccccccccc

      open (unit=10,file=
     &     'CAGE.pdb',
     &     status='old')

      do i=1,num_SubUnit_template
         do j=1,npart
            read(10,1200) res_x(j,i),res_y(j,i),res_z(j,i)
         enddo
c         read(10,*)
      enddo

      close(10)

 1200 format(30x,3F8.3)
 

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   generate multiple trajectories
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do n_t=1,num_trajec

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c   construct the initialized position and conformation of molecules in 3D
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   construct the position and conformation of each molecule

c>>>>>   random position

         do i=1,mol_tot_num
 100        continue
            temp_i=rand3(r3)*cell_range_x-cell_range_x/2  
            temp_j=rand3(r3)*cell_range_y-cell_range_y/2  
            temp_k=rand3(r3)*cell_range_z-cell_range_z/2            
            do j=1,i-1
               dist=sqrt((temp_i-mol_x(npart,j))**2+
     &              (temp_j-mol_y(npart,j))**2+
     &              (temp_k-mol_z(npart,j))**2)
               if(dist.le.
     &              2*radius)then
                  goto 100
               endif
            enddo
            do j=1,npart
               mol_x_0(j,i)=res_x(j,1)-(res_x(6,1)
     &              +res_x(7,1))*0.5+temp_i
               mol_y_0(j,i)=res_y(j,1)-(res_y(6,1)
     &              +res_y(7,1))*0.5+temp_j
               mol_z_0(j,i)=res_z(j,1)-(res_z(6,1)
     &              +res_z(7,1))*0.5+temp_k
            enddo
            mol_status(i)=0
            do j=1,npart
               site_status(j,i)=0
            enddo
            mol_type_index(i)=int(rand3(r3)*mol_type_num)+1
            
c>>>>>   random orientation

            theta=(2*rand3(r3)-1)*pai
            phi=(2*rand3(r3)-1)*pai
            psai=(2*rand3(r3)-1)*pai
            
            t(1,1)=cos(psai)*cos(phi)-cos(theta)*sin(phi)*sin(psai)
            t(1,2)=-sin(psai)*cos(phi)-cos(theta)*sin(phi)*cos(psai)
            t(1,3)=sin(theta)*sin(phi)
            
            t(2,1)=cos(psai)*sin(phi)+cos(theta)*cos(phi)*sin(psai)
            t(2,2)=-sin(psai)*sin(phi)+cos(theta)*cos(phi)*cos(psai)
            t(2,3)=-sin(theta)*cos(phi)
            
            t(3,1)=sin(psai)*sin(theta)
            t(3,2)=cos(psai)*sin(theta)
            t(3,3)=cos(theta)
         
            do k=1,npart
               mol_x0(k,i)=t(1,1)*(mol_x_0(k,i)-
     &              (mol_x_0(6,i)+mol_x_0(7,i))*0.5)
     &              +t(1,2)*(mol_y_0(k,i)-(mol_y_0(6,i)
     &              +mol_y_0(7,i))*0.5)
     &              +t(1,3)*(mol_z_0(k,i)-(mol_z_0(6,i)
     &              +mol_z_0(7,i))*0.5)
     &              +(mol_x_0(6,i)+mol_x_0(7,i))*0.5
               mol_y0(k,i)=t(2,1)*(mol_x_0(k,i)-
     &              (mol_x_0(6,i)+mol_x_0(7,i))*0.5)
     &              +t(2,2)*(mol_y_0(k,i)-(mol_y_0(6,i)
     &              +mol_y_0(7,i))*0.5)
     &              +t(2,3)*(mol_z_0(k,i)-(mol_z_0(6,i)
     &              +mol_z_0(7,i))*0.5)
     &              +(mol_y_0(6,i)+mol_y_0(7,i))*0.5
               mol_z0(k,i)=t(3,1)*(mol_x_0(k,i)-
     &              (mol_x_0(6,i)+mol_x_0(7,i))*0.5)
     &              +t(3,2)*(mol_y_0(k,i)-(mol_y_0(6,i)
     &              +mol_y_0(7,i))*0.5)
     &              +t(3,3)*(mol_z_0(k,i)-(mol_z_0(6,i)
     &              +mol_z_0(7,i))*0.5)
     &              +(mol_z_0(6,i)+mol_z_0(7,i))*0.5
            enddo
         
         enddo

         do i=1,mol_tot_num
            do j=1,npart
               mol_x(j,i)=mol_x0(j,i)
               mol_y(j,i)=mol_y0(j,i)
               mol_z(j,i)=mol_z0(j,i)		 
            enddo
         enddo
         
         complex_num=0
         do i=1,mol_tot_num*npart
            do j=1,2
               complex_mol_idx(i,j)=0
               complex_site_idx(i,j)=0
            enddo
         enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      begin  main  loop of Diffusion-Reaction simulation
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         initial_simu_time=0.0

         do mc_time_step=1,simu_step

            current_simu_time=current_simu_time+time_step
c            print*,mc_time_step
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
            do i=1,mol_tot_num
               do j=1,npart
                  mol_x_new(j,i)=mol_x(j,i)
                  mol_y_new(j,i)=mol_y(j,i)
                  mol_z_new(j,i)=mol_z(j,i)
                  mol_x0_new(j,i)=mol_x0(j,i)
                  mol_y0_new(j,i)=mol_y0(j,i)
                  mol_z0_new(j,i)=mol_z0(j,i)				  
                  site_status_new(j,i)=site_status(j,i)
               enddo
               mol_status_new(i)=0
               do j=1,npart
                  if(site_status_new(j,i).eq.1)then
                     mol_status_new(i)=1
                  endif
               enddo
            enddo
            complex_num_new=complex_num
            do i=1,complex_num
               do j=1,2
                  complex_mol_idx_new(i,j)=complex_mol_idx(i,j)
                  complex_site_idx_new(i,j)=complex_site_idx(i,j)
               enddo
            enddo

            do iteration_mole_step=1,mol_tot_num

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   randomly select one molecule
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               mol_index=int(rand3(r3)*mol_tot_num)+1

cc>>>>   make random diffusion for selected molecule
c>>>>   if it is unbound
               if(mol_status_new(mol_index).eq.0)then
                  Prob_Diff=(6*mol_D*time_step)/(distance_step**2)
                  if(rand3(r3).lt.Prob_Diff)then
cc>>     if smaller than the diffusion probability, move alone
                     theta=rand3(r3)*pai
                     phai=rand3(r3)*2*pai
                     displace_amp=2*distance_step*rand3(r3)
                     do i=1,npart
                        mol_x0_new0(i,mol_index)=
     &                       mol_x0(i,mol_index)+
     &                       displace_amp*sin(theta)*cos(phai)
                        mol_y0_new0(i,mol_index)=
     &                       mol_y0(i,mol_index)+
     &                       displace_amp*sin(theta)*sin(phai)
                        mol_z0_new0(i,mol_index)=
     &                       mol_z0(i,mol_index)+
     &                       displace_amp*cos(theta)
                     enddo
                     PB_x=cell_range_x*
     &                    anint((mol_x0_new0(6,mol_index)+
     &                    mol_x0_new0(7,mol_index))*0.5/cell_range_x)
                     PB_y=cell_range_y*
     &                    anint((mol_y0_new0(6,mol_index)+
     &                    mol_y0_new0(7,mol_index))*0.5/cell_range_y)
                     PB_z=cell_range_z*
     &                    anint((mol_z0_new0(6,mol_index)+
     &                    mol_z0_new0(7,mol_index))*0.5/cell_range_z)
                     do i=1,npart
                        mol_x0_new0(i,mol_index)=
     &                       mol_x0_new0(i,mol_index)-PB_x     
                        mol_y0_new0(i,mol_index)=
     &                       mol_y0_new0(i,mol_index)-PB_y
                        mol_z0_new0(i,mol_index)=
     &                       mol_z0_new0(i,mol_index)-PB_z
                     enddo
                  else
                     do i=1,npart
                        mol_x0_new0(i,mol_index)=
     &                       mol_x0(i,mol_index)
                        mol_y0_new0(i,mol_index)=
     &                       mol_y0(i,mol_index)
                        mol_z0_new0(i,mol_index)=
     &                       mol_z0(i,mol_index)
                     enddo
                  endif

c>>>>>>>>> Rotate the molecule

                  theta=(2*rand3(r3)-1)*mol_rot_D*pai/180.0
                  phi=(2*rand3(r3)-1)*mol_rot_D*pai/180.0 
                  psai=(2*rand3(r3)-1)*mol_rot_D*pai/180.0
                  
                  t(1,1)=cos(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*sin(psai)
                  t(1,2)=-sin(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*cos(psai)
                  t(1,3)=sin(theta)*sin(phi)
                  
                  t(2,1)=cos(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*sin(psai)
                  t(2,2)=-sin(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*cos(psai)
                  t(2,3)=-sin(theta)*cos(phi)
                  
                  t(3,1)=sin(psai)*sin(theta)
                  t(3,2)=cos(psai)*sin(theta)
                  t(3,3)=cos(theta)

                  
                  do k=1,npart
                     mol_x0_new(k,mol_index)=t(1,1)*
     &                    (mol_x0_new0(k,mol_index)
     &                    -(mol_x0_new0(6,mol_index)+
     &                    mol_x0_new0(7,mol_index))*0.5)+
     &                    t(1,2)*(mol_y0_new0(k,mol_index)-
     &                    (mol_y0_new0(6,mol_index)+
     &                    mol_y0_new0(7,mol_index))*0.5)
     &                    +t(1,3)*(mol_z0_new0(k,mol_index)-
     &                    (mol_z0_new0(6,mol_index)+
     &                    mol_z0_new0(7,mol_index))*0.5)
     &                    +(mol_x0_new0(6,mol_index)+
     &                    mol_x0_new0(7,mol_index))*0.5
                     mol_y0_new(k,mol_index)=t(2,1)*
     &                    (mol_x0_new0(k,mol_index)
     &                    -(mol_x0_new0(6,mol_index)+
     &                    mol_x0_new0(7,mol_index))*0.5)+
     &                    t(2,2)*(mol_y0_new0(k,mol_index)-
     &                    (mol_y0_new0(6,mol_index)+
     &                    mol_y0_new0(7,mol_index))*0.5)
     &                    +t(2,3)*(mol_z0_new0(k,mol_index)-
     &                    (mol_z0_new0(6,mol_index)+
     &                    mol_z0_new0(7,mol_index))*0.5)
     &                    +(mol_y0_new0(6,mol_index)+
     &                    mol_y0_new0(7,mol_index))*0.5
                     mol_z0_new(k,mol_index)=t(3,1)*
     &                    (mol_x0_new0(k,mol_index)
     &                    -(mol_x0_new0(6,mol_index)+
     &                    mol_x0_new0(7,mol_index))*0.5)+
     &                    t(3,2)*(mol_y0_new0(k,mol_index)-
     &                    (mol_y0_new0(6,mol_index)+
     &                    mol_y0_new0(7,mol_index))*0.5)
     &                    +t(3,3)*(mol_z0_new0(k,mol_index)-
     &                    (mol_z0_new0(6,mol_index)+
     &                    mol_z0_new0(7,mol_index))*0.5)
     &                    +(mol_z0_new0(6,mol_index)+
     &                    mol_z0_new0(7,mol_index))*0.5
                  enddo
				  
c>>>>>   random translational flexibility


                  Flx_Rg_x=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  Flx_Rg_y=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  Flx_Rg_z=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  mol_x_new0(1,mol_index)=
     &                 mol_x0_new(1,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(1,mol_index)=
     &                 mol_y0_new(1,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(1,mol_index)=
     &                 mol_z0_new(1,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(3,mol_index)=
     &                 mol_x0_new(3,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(3,mol_index)=
     &                 mol_y0_new(3,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(3,mol_index)=
     &                 mol_z0_new(3,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(4,mol_index)=
     &                 mol_x0_new(4,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(4,mol_index)=
     &                 mol_y0_new(4,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(4,mol_index)=
     &                 mol_z0_new(4,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(6,mol_index)=
     &                 mol_x0_new(6,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(6,mol_index)=
     &                 mol_y0_new(6,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(6,mol_index)=
     &                 mol_z0_new(6,mol_index)
     &                 +Flx_Rg_z
                  
                  Flx_Rg_x=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  Flx_Rg_y=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  Flx_Rg_z=2*Flx_TlRg*rand4(r4)-Flx_TlRg
                  mol_x_new0(2,mol_index)=
     &                 mol_x0_new(2,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(2,mol_index)=
     &                 mol_y0_new(2,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(2,mol_index)=
     &                 mol_z0_new(2,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(5,mol_index)=
     &                 mol_x0_new(5,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(5,mol_index)=
     &                 mol_y0_new(5,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(5,mol_index)=
     &                 mol_z0_new(5,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(7,mol_index)=
     &                 mol_x0_new(7,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(7,mol_index)=
     &                 mol_y0_new(7,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(7,mol_index)=
     &                 mol_z0_new(7,mol_index)
     &                 +Flx_Rg_z
                  mol_x_new0(8,mol_index)=
     &                 mol_x0_new(8,mol_index)
     &                 +Flx_Rg_x
                  mol_y_new0(8,mol_index)=
     &                 mol_y0_new(8,mol_index)
     &                 +Flx_Rg_y 
                  mol_z_new0(8,mol_index)=
     &                 mol_z0_new(8,mol_index)
     &                 +Flx_Rg_z	 
                  
c>>>>>   random rotational flexibility

                  theta=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  phi=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  psai=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  
                  t(1,1)=cos(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*sin(psai)
                  t(1,2)=-sin(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*cos(psai)
                  t(1,3)=sin(theta)*sin(phi)
                  
                  t(2,1)=cos(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*sin(psai)
                  t(2,2)=-sin(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*cos(psai)
                  t(2,3)=-sin(theta)*cos(phi)
                  
                  t(3,1)=sin(psai)*sin(theta)
                  t(3,2)=cos(psai)*sin(theta)
                  t(3,3)=cos(theta)
                  
                  do k=1,npart
                     if((k.eq.1).or.(k.eq.3).or.
     &                    (k.eq.4).or.(k.eq.6))then
                        mol_x_new(k,mol_index)=t(1,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(1,mol_index)+
     &                       mol_x_new0(1,mol_index))*0.5)+
     &                       t(1,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(1,mol_index)+
     &                       mol_y_new0(1,mol_index))*0.5)
     &                       +t(1,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(1,mol_index)+
     &                       mol_z_new0(1,mol_index))*0.5)
     &                       +(mol_x_new0(1,mol_index)+
     &                       mol_x_new0(1,mol_index))*0.5
                        mol_y_new(k,mol_index)=t(2,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(1,mol_index)+
     &                       mol_x_new0(1,mol_index))*0.5)+
     &                       t(2,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(1,mol_index)+
     &                       mol_y_new0(1,mol_index))*0.5)
     &                       +t(2,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(1,mol_index)+
     &                       mol_z_new0(1,mol_index))*0.5)
     &                       +(mol_y_new0(1,mol_index)+
     &                       mol_y_new0(1,mol_index))*0.5
                        mol_z_new(k,mol_index)=t(3,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(1,mol_index)+
     &                       mol_x_new0(1,mol_index))*0.5)+
     &                       t(3,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(1,mol_index)+
     &                       mol_y_new0(1,mol_index))*0.5)
     &                       +t(3,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(1,mol_index)+
     &                       mol_z_new0(1,mol_index))*0.5)
     &                       +(mol_z_new0(1,mol_index)+
     &                       mol_z_new0(1,mol_index))*0.5	 
                     endif 
                  enddo
                  
                  theta=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  phi=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  psai=(2*rand5(r5)-1)*Flx_RtRg*pai/180.0
                  
                  t(1,1)=cos(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*sin(psai)
                  t(1,2)=-sin(psai)*cos(phi)-
     &                 cos(theta)*sin(phi)*cos(psai)
                  t(1,3)=sin(theta)*sin(phi)
                  
                  t(2,1)=cos(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*sin(psai)
                  t(2,2)=-sin(psai)*sin(phi)+
     &                 cos(theta)*cos(phi)*cos(psai)
                  t(2,3)=-sin(theta)*cos(phi)
                  
                  t(3,1)=sin(psai)*sin(theta)
                  t(3,2)=cos(psai)*sin(theta)
                  t(3,3)=cos(theta)
                  
                  do k=1,npart
                     if((k.eq.2).or.(k.eq.5).or.
     &                    (k.eq.7).or.(k.eq.8))then
                        mol_x_new(k,mol_index)=t(1,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(2,mol_index)+
     &                       mol_x_new0(2,mol_index))*0.5)+
     &                       t(1,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(2,mol_index)+
     &                       mol_y_new0(2,mol_index))*0.5)
     &                       +t(1,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(2,mol_index)+
     &                       mol_z_new0(2,mol_index))*0.5)
     &                       +(mol_x_new0(2,mol_index)+
     &                       mol_x_new0(2,mol_index))*0.5
                        mol_y_new(k,mol_index)=t(2,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(2,mol_index)+
     &                       mol_x_new0(2,mol_index))*0.5)+
     &                       t(2,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(2,mol_index)+
     &                       mol_y_new0(2,mol_index))*0.5)
     &                       +t(2,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(2,mol_index)+
     &                       mol_z_new0(2,mol_index))*0.5)
     &                       +(mol_y_new0(2,mol_index)+
     &                       mol_y_new0(2,mol_index))*0.5
                        mol_z_new(k,mol_index)=t(3,1)*
     &                       (mol_x_new0(k,mol_index)
     &                       -(mol_x_new0(2,mol_index)+
     &                       mol_x_new0(2,mol_index))*0.5)+
     &                       t(3,2)*(mol_y_new0(k,mol_index)-
     &                       (mol_y_new0(2,mol_index)+
     &                       mol_y_new0(2,mol_index))*0.5)
     &                       +t(3,3)*(mol_z_new0(k,mol_index)-
     &                       (mol_z_new0(2,mol_index)+
     &                       mol_z_new0(2,mol_index))*0.5)
     &                       +(mol_z_new0(2,mol_index)+
     &                       mol_z_new0(2,mol_index))*0.5	 
                     endif 
                  enddo
						
						
ccc>>>>   Check Collisions

                  collision_flag=0
                  do i=1,mol_tot_num
                     if(mol_index.ne.i)then
                        dist=sqrt(((mol_x_new(6,i)+mol_x_new(7,i))*0.5-
     &                       (mol_x_new(6,mol_index)+
     &                       mol_x_new(7,mol_index))*0.5)**2+
     &                       ((mol_y_new(6,i)+mol_y_new(7,i))*0.5-
     &                       (mol_y_new(6,mol_index)+
     &                       mol_y_new(7,mol_index))*0.5)**2+
     &                       ((mol_z_new(6,i)+mol_z_new(7,i))*0.5-
     &                       (mol_z_new(6,mol_index)+
     &                       mol_z_new(7,mol_index))*0.5)**2)
                        if(dist.lt.radius*2)then
                           collision_flag=1
                        endif
                     endif
                  enddo
                  
                  if(collision_flag.eq.1)then
                     do i=1,npart
                        mol_x_new0(i,mol_index)=
     &                       mol_x(i,mol_index)
                        mol_y_new0(i,mol_index)=
     &                       mol_y(i,mol_index)
                        mol_z_new0(i,mol_index)=
     &                       mol_z(i,mol_index)
                        mol_x_new(i,mol_index)=
     &                       mol_x(i,mol_index)
                        mol_y_new(i,mol_index)=
     &                       mol_y(i,mol_index)
                        mol_z_new(i,mol_index)=
     &                       mol_z(i,mol_index)
                        mol_x0_new0(i,mol_index)=
     &                       mol_x0(i,mol_index)
                        mol_y0_new0(i,mol_index)=
     &                       mol_y0(i,mol_index)
                        mol_z0_new0(i,mol_index)=
     &                       mol_z0(i,mol_index)
                        mol_x0_new(i,mol_index)=
     &                       mol_x0(i,mol_index)
                        mol_y0_new(i,mol_index)=
     &                       mol_y0(i,mol_index)
                        mol_z0_new(i,mol_index)=
     &                       mol_z0(i,mol_index)	 
                     enddo
                  endif

               endif

            enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   endo of diffusion for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
c            goto 1977
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   perform the reaction for molecules
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccc  1) Two subunits associate into a dimer

            do i=1,mol_tot_num-1
               do j=i+1,mol_tot_num

                  do i2=1,npart
                     do j2=1,npart

                        if((site_status_new(i2,i).eq.0)
     &                       .AND.(site_status_new(j2,j).eq.0))then
                           
c>>>>>>>>>>>>>   calculate distance of two reaction sites between two subunits

                           dist=sqrt((mol_x_new(j2,j)-
     &                          mol_x_new(i2,i))**2+
     &                          (mol_y_new(j2,j)-
     &                          mol_y_new(i2,i))**2+
     &                          (mol_z_new(j2,j)-
     &                          mol_z_new(i2,i))**2)
                           
                           if((dist.lt.bond_dist_cutoff).AND.
     &                          (((i2.eq.3).AND.(j2.eq.4)).OR.
     &                          ((i2.eq.4).AND.(j2.eq.3)).OR.
     &                          ((i2.eq.5).AND.(j2.eq.5))))then
                              
                              
                              if(((i2.eq.3).AND.
     &                             (j2.eq.4)).OR.
     &                             ((i2.eq.4).AND.(j2.eq.3)))then
	 
                         
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                                 theta_pd=0
                                 point_x(1)=mol_x_new(i2,i)
     &                                -mol_x_new(1,i)
                                 point_y(1)=mol_y_new(i2,i)
     &                                -mol_y_new(1,i)
                                 point_z(1)=mol_z_new(i2,i)
     &                                -mol_z_new(1,i)
                                 point_x(2)=0
                                 point_y(2)=0
                                 point_z(2)=0
                                 point_x(3)=mol_x_new(j2,j)
     &                                -mol_x_new(1,j)
                                 point_y(3)=mol_y_new(j2,j)
     &                                -mol_y_new(1,j)
                                 point_z(3)=mol_z_new(j2,j)
     &                                -mol_z_new(1,j)
                                 call gettheta(point_x,point_y,point_z,
     &                                theta_pd)
c>> orientational angle
                                 theta_ot=0
                                 point_x(1)=mol_x_new(6,i)
     &                                -mol_x_new(1,i)
                                 point_y(1)=mol_y_new(6,i)
     &                                -mol_y_new(1,i)
                                 point_z(1)=mol_z_new(6,i)
     &                                -mol_z_new(1,i)
                                 point_x(2)=0
                                 point_y(2)=0
                                 point_z(2)=0
                                 point_x(3)=mol_x_new(6,j)
     &                                -mol_x_new(1,j)
                                 point_y(3)=mol_y_new(6,j)
     &                                -mol_y_new(1,j)
                                 point_z(3)=mol_z_new(6,j)
     &                                -mol_z_new(1,j)
                                 call gettheta(point_x,point_y,point_z,
     &                                theta_ot)
c                           print*,theta_pd,theta_ot
                                 if((abs(theta_ot-bond_thetaot_1)
     &                                .lt.bond_thetaot_cutoff_1).AND.
     &                                (abs(theta_pd-bond_thetapd_1)
     &                                .lt.bond_thetapd_cutoff_1))then
                                    
                                    Ass_Rate_Value=Ass_Rate_max
     &                                   (mol_type_index(i),
     &                                   mol_type_index(j),
     &                                   1)
                              
                                    Prob_Ass=Ass_Rate_Value*time_step ! ready to change
                                    if(rand3(r3).lt.Prob_Ass)then	 
	 
                                
                                       if((mol_status_new(i).eq.1).AND.
     &                                      (mol_status_new(j).eq.1)
     &                                      )then
                                									   
                                       
                                          site_status_new(i2,i)=1
                                          site_status_new(j2,j)=1
                                          complex_num_new=
     &                                         complex_num_new+1
                                          complex_mol_idx_new
     &                                         (complex_num_new,1)=i
                                          complex_mol_idx_new
     &                                         (complex_num_new,2)=j
                                          complex_site_idx_new
     &                                         (complex_num_new,1)=i2
                                          complex_site_idx_new
     &                                         (complex_num_new,2)=j2
                                          
                                       else
                                    
                                          align_num_i=0
                                          align_num_j=0
                                          target_ind_i2=0
                                          target_ind_j2=0
                                          
                                          if((mol_status_new(i).eq.0)
     &                                         .AND.
     &                                         (mol_status_new(j).eq.0))
     &                                         then
                                             align_num_i=npart
                                             align_num_j=npart
                                             align_ind=0
                                          elseif(
     &                                            (mol_status_new(i)
     &                                            .eq.1)
     &                                            .AND.
     &                                            (mol_status_new(j)
     &                                            .eq.0))
     &                                            then
                                             align_num_i=npart
                                             align_num_j=0
                                             align_ind=0
                                          elseif((mol_status_new(i)
     &                                            .eq.0)
     &                                            .AND.
     &                                            (mol_status_new(j)
     &                                            .eq.1))
     &                                            then
                                             align_num_i=0
                                             align_num_j=npart
                                             align_ind=npart
                                          endif
                                    
                                          if((i2.eq.3).AND.(j2.eq.4)
     &                                         )then
                                             target_ind_i2=8
                                             target_ind_j2=10
                                          elseif((i2.eq.4).AND.
     &                                            (j2.eq.3))then
                                             target_ind_i2=8
                                             target_ind_j2=12
                                          endif
                                    
                                          do k=1,1000
                                             xa(k)=0
                                             ya(k)=0
                                             za(k)=0
                                             xb(k)=0
                                             yb(k)=0
                                             zb(k)=0
                                             atompr(1,k)=0
                                             atompr(2,k)=0
                                          enddo
                                          index=0
                                          do k=1,align_num_i
                                             index=index+1
                                             xa(index)=mol_x_new(k,i)
                                             ya(index)=mol_y_new(k,i)
                                             za(index)=mol_z_new(k,i)
                                          enddo
                                          do k=1,align_num_j
                                             index=index+1
                                             xa(index)=mol_x_new(k,j)
                                             ya(index)=mol_y_new(k,j)
                                             za(index)=mol_z_new(k,j)
                                          enddo
                                          index=0
                                          do k=1,npart
                                             index=index+1
                                             xb(index)=res_x
     &                                            (k,target_ind_i2)
                                             yb(index)=res_y
     &                                            (k,target_ind_i2)
                                             zb(index)=res_z
     &                                            (k,target_ind_i2)
                                          enddo
                                          do k=1,npart
                                             index=index+1
                                             xb(index)=res_x
     &                                            (k,target_ind_j2)
                                             yb(index)=res_y
     &                                            (k,target_ind_j2)
                                             zb(index)=res_z
     &                                            (k,target_ind_j2)
                                          enddo
                                          do k=1,align_num_i+align_num_j
                                             atompr(1,k)=k
                                             atompr(2,k)=k+align_ind
                                          enddo
                                          npair=align_num_i+align_num_j       
                                          npair_a=align_num_i+
     &                                         align_num_j
                                          npair_b=npart*2
                                          call ROTLSQ
     &                                         (xa,ya,za,npair_a,
     &                                         xb,yb,zb,
     &                                         npair_b,atompr,npair)
                                    
c>>  check collision

                                          collision_flag=0
                                          do k=1,mol_tot_num
                                             if((k.ne.i).AND.(k.ne.j)
     &                                            )then
                                                dist=sqrt(((mol_x_new
     &                                               (6,k)+
     &                                               mol_x_new(7,k))
     &                                               *0.5-
     &                                               (xb(6)+xb(7))
     &                                               *0.5)**2+
     &                                               ((mol_y_new
     &                                               (6,k)+
     &                                               mol_y_new(7,k))
     &                                               *0.5-
     &                                               (yb(6)+yb(7))
     &                                               *0.5)**2+
     &                                               ((mol_z_new
     &                                               (6,k)+
     &                                               mol_z_new(7,k))
     &                                               *0.5-
     &                                               (zb(6)+zb(7))
     &                                               *0.5)**2)
                                                if(dist.lt.radius*2)then
                                                   collision_flag=1
                                                endif
                                             endif
                                          enddo
 
                                          do k=1,mol_tot_num
                                             if((k.ne.i).AND.(k.ne.j)
     &                                            )then
                                                dist=sqrt(((mol_x_new
     &                                               (6,k)+
     &                                               mol_x_new(7,k))
     &                                               *0.5-
     &                                               (xb(14)+xb(15))
     &                                               *0.5)**2+
     &                                               ((mol_y_new
     &                                               (6,k)+
     &                                               mol_y_new(7,k))
     &                                               *0.5-
     &                                               (yb(14)+yb(15))
     &                                               *0.5)**2+
     &                                               ((mol_z_new
     &                                               (6,k)+
     &                                               mol_z_new(7,k))
     &                                               *0.5-
     &                                               (zb(14)+zb(15))
     &                                               *0.5)**2)
                                                if(dist.lt.radius*2)then
                                                   collision_flag=1
                                                endif
                                             endif
                                          enddo
  
c>>   if no collision, reaction
                                    
                                          if(collision_flag.eq.0)then
                                             
                                             site_status_new(i2,i)=1
                                             site_status_new(j2,j)=1
                                             complex_num_new=
     &                                            complex_num_new+1
                                             complex_mol_idx_new
     &                                      (complex_num_new,1)=i
                                             complex_mol_idx_new
     &                                            (complex_num_new,2)=j
                                             complex_site_idx_new
     &                                            (complex_num_new,1)=i2
                                             complex_site_idx_new
     &                                            (complex_num_new,2)=j2
                                             
                                             index=0
                                             do k=1,npart
                                                index=index+1
                                                mol_x_new(k,i)=xb(index)
                                                mol_y_new(k,i)=yb(index)                                                
                                                mol_z_new(k,i)=zb(index)
                                                mol_x0_new(k,i)=
     &                                               xb(index)
                                                mol_y0_new(k,i)=
     &                                               yb(index)                                                
                                                mol_z0_new(k,i)=
     &                                               zb(index)												
                                             enddo
                                             do k=1,npart
                                                index=index+1
                                                mol_x_new(k,j)=xb(index)
                                                mol_y_new(k,j)=yb(index)
                                                mol_z_new(k,j)=zb(index)
                                                mol_x0_new(k,j)=
     &                                               xb(index)
                                                mol_y0_new(k,j)=
     &                                               yb(index)
                                                mol_z0_new(k,j)=
     &                                               zb(index)												
                                             enddo
                                             
                                             
                                             
c>>   after reaction update mol_status
                                       
                                             do k=1,mol_tot_num
                                                mol_status_new(k)=0
                                                do k2=1,npart
                                                   if(site_status_new
     &                                                  (k2,k)
     &                                                  .eq.1)then
                                                      mol_status_new(k)
     &                                                     =1
                                                   endif
                                                enddo
                                             enddo
                                             

                                          endif                                       
                                          
                                       endif	   
                                    endif

                                 endif                                       
                                   
                                 
                              elseif((i2.eq.5).AND.
     &                                (j2.eq.5))then
                                 
c>>>>>>>>>>>> calculate two angle
c>> perpendicular angle
                                 theta_pd=0
                                 point_x(1)=mol_x_new(i2,i)
     &                                -mol_x_new(2,i)
                                 point_y(1)=mol_y_new(i2,i)
     &                                -mol_y_new(2,i)
                                 point_z(1)=mol_z_new(i2,i)
     &                                -mol_z_new(2,i)
                                 point_x(2)=0
                                 point_y(2)=0
                                 point_z(2)=0
                                 point_x(3)=mol_x_new(j2,j)
     &                                -mol_x_new(2,j)
                                 point_y(3)=mol_y_new(j2,j)
     &                                -mol_y_new(2,j)
                                 point_z(3)=mol_z_new(j2,j)
     &                                -mol_z_new(2,j)
                                 call gettheta(point_x,point_y,point_z,
     &                                theta_pd)
c>>   orientational angle
                                 theta_ot=0
                                 point_x(1)=mol_x_new(8,i)
     &                                -mol_x_new(2,i)
                                 point_y(1)=mol_y_new(8,i)
     &                                -mol_y_new(2,i)
                                 point_z(1)=mol_z_new(8,i)
     &                                -mol_z_new(2,i)
                                 point_x(2)=0
                                 point_y(2)=0
                                 point_z(2)=0
                                 point_x(3)=mol_x_new(8,j)
     &                                -mol_x_new(2,j)
                                 point_y(3)=mol_y_new(8,j)
     &                                -mol_y_new(2,j)
                                 point_z(3)=mol_z_new(8,j)
     &                                -mol_z_new(2,j)
                                 call gettheta(point_x,point_y,point_z,
     &                                theta_ot)
c                           print*,theta_pd,theta_ot
                                 if((abs(theta_ot-bond_thetaot_2)
     &                                .lt.bond_thetaot_cutoff_2).AND.
     &                                (abs(theta_pd-bond_thetapd_2)
     &                                .lt.bond_thetapd_cutoff_2))then
                                    
                                    Ass_Rate_Value=Ass_Rate_max
     &                                   (mol_type_index(i),
     &                                   mol_type_index(j),
     &                                   2)
                                    
                                    Prob_Ass=Ass_Rate_Value*time_step ! ready to change
                                    if(rand3(r3).lt.Prob_Ass)then	 
                                       
                                 
                                       if((mol_status_new(i).eq.1).AND.
     &                                      (mol_status_new(j).eq.1)
     &                                      )then
                                          
                                       
                                          site_status_new(i2,i)=1
                                          site_status_new(j2,j)=1
                                          complex_num_new=
     &                                         complex_num_new+1
                                          complex_mol_idx_new
     &                                         (complex_num_new,1)=i
                                          complex_mol_idx_new
     &                                         (complex_num_new,2)=j
                                          complex_site_idx_new
     &                                         (complex_num_new,1)=i2
                                          complex_site_idx_new
     &                                         (complex_num_new,2)=j2
                                          
                                       else
                                          
                                          align_num_i=0
                                          align_num_j=0
                                          target_ind_i2=0
                                          target_ind_j2=0
                                          
                                          if((mol_status_new(i).eq.0)
     &                                         .AND.
     &                                         (mol_status_new(j).eq.0))
     &                                         then
                                             align_num_i=npart
                                             align_num_j=npart
                                             align_ind=0
                                          elseif((mol_status_new(i)
     &                                            .eq.1)
     &                                            .AND.
     &                                            (mol_status_new(j)
     &                                            .eq.0))
     &                                            then
                                             align_num_i=npart
                                             align_num_j=0
                                             align_ind=0
                                          elseif((mol_status_new(i)
     &                                            .eq.0)
     &                                            .AND.
     &                                            (mol_status_new(j)
     &                                            .eq.1))
     &                                            then
                                             align_num_i=0
                                             align_num_j=npart
                                             align_ind=npart
                                          endif
                                          
                                          if((i2.eq.5).AND.(j2.eq.5)
     &                                         )then
                                             target_ind_i2=8
                                             target_ind_j2=7
                                          endif
                                          
                                          do k=1,1000
                                             xa(k)=0
                                             ya(k)=0
                                             za(k)=0
                                             xb(k)=0
                                             yb(k)=0
                                             zb(k)=0
                                             atompr(1,k)=0
                                             atompr(2,k)=0
                                          enddo
                                          index=0
                                          do k=1,align_num_i
                                             index=index+1
                                             xa(index)=mol_x_new(k,i)
                                             ya(index)=mol_y_new(k,i)
                                             za(index)=mol_z_new(k,i)
                                          enddo
                                          do k=1,align_num_j
                                             index=index+1
                                             xa(index)=mol_x_new(k,j)
                                             ya(index)=mol_y_new(k,j)
                                             za(index)=mol_z_new(k,j)
                                          enddo
                                          index=0
                                          do k=1,npart
                                             index=index+1
                                             xb(index)=res_x
     &                                            (k,target_ind_i2)
                                             yb(index)=res_y
     &                                            (k,target_ind_i2)
                                             zb(index)=res_z
     &                                            (k,target_ind_i2)
                                          enddo
                                          do k=1,npart
                                             index=index+1
                                             xb(index)=res_x
     &                                            (k,target_ind_j2)
                                             yb(index)=res_y
     &                                            (k,target_ind_j2)
                                             zb(index)=res_z
     &                                            (k,target_ind_j2)
                                          enddo
                                          do k=1,align_num_i+align_num_j
                                             atompr(1,k)=k
                                             atompr(2,k)=k+align_ind
                                          enddo
                                          npair=align_num_i+align_num_j       
                                          npair_a=align_num_i+
     &                                         align_num_j
                                          npair_b=npart*2
                                          call ROTLSQ
     &                                         (xa,ya,za,npair_a,
     &                                         xb,yb,zb,
     &                                         npair_b,atompr,npair)
                                          
c>>  check collision

                                          collision_flag=0
                                          do k=1,mol_tot_num
                                             if((k.ne.i).AND.(k.ne.j)
     &                                            )then
                                                dist=sqrt(((mol_x_new
     &                                               (6,k)+
     &                                               mol_x_new(7,k))
     &                                               *0.5-
     &                                               (xb(6)+xb(7))
     &                                               *0.5)**2+
     &                                               ((mol_y_new
     &                                               (6,k)+
     &                                               mol_y_new(7,k))
     &                                               *0.5-
     &                                               (yb(6)+yb(7))
     &                                               *0.5)**2+
     &                                               ((mol_z_new
     &                                               (6,k)+
     &                                               mol_z_new(7,k))
     &                                               *0.5-
     &                                               (zb(6)+zb(7))
     &                                               *0.5)**2)
                                                if(dist.lt.radius*2)then
                                                   collision_flag=1
                                                endif
                                             endif
                                          enddo
 
                                          do k=1,mol_tot_num
                                             if((k.ne.i).AND.(k.ne.j)
     &                                            )then
                                                dist=sqrt(((mol_x_new
     &                                               (6,k)+
     &                                               mol_x_new(7,k))
     &                                               *0.5-
     &                                               (xb(14)+xb(15))
     &                                               *0.5)**2+
     &                                               ((mol_y_new
     &                                               (6,k)+
     &                                               mol_y_new(7,k))
     &                                               *0.5-
     &                                               (yb(14)+yb(15))
     &                                               *0.5)**2+
     &                                               ((mol_z_new
     &                                               (6,k)+
     &                                               mol_z_new(7,k))
     &                                               *0.5-
     &                                               (zb(14)+zb(15))
     &                                               *0.5)**2)
                                                if(dist.lt.radius*2)then
                                                   collision_flag=1
                                                endif
                                             endif
                                          enddo
              
c>>   if no collision, reaction
                                          
                                          if(collision_flag.eq.0)then
                                             
                                             site_status_new(i2,i)=1
                                             site_status_new(j2,j)=1
                                             complex_num_new=
     &                                            complex_num_new+1
                                             complex_mol_idx_new
     &                                            (complex_num_new,1)=i
                                             complex_mol_idx_new
     &                                            (complex_num_new,2)=j
                                             complex_site_idx_new
     &                                            (complex_num_new,1)=i2
                                             complex_site_idx_new
     &                                            (complex_num_new,2)=j2
                                             
                                             index=0
                                             do k=1,npart
                                                index=index+1
                                                mol_x_new(k,i)=xb(index)
                                                mol_y_new(k,i)=yb(index)                                                
                                                mol_z_new(k,i)=zb(index)
                                                mol_x0_new(k,i)=
     &                                               xb(index)
                                                mol_y0_new(k,i)=
     &                                               yb(index)                                                
                                                mol_z0_new(k,i)=
     &                                               zb(index)												
                                             enddo
                                             do k=1,npart
                                                index=index+1
                                                mol_x_new(k,j)=xb(index)
                                                mol_y_new(k,j)=yb(index)
                                                mol_z_new(k,j)=zb(index)
                                                mol_x0_new(k,j)=
     &                                               xb(index)
                                                mol_y0_new(k,j)=
     &                                               yb(index)
                                                mol_z0_new(k,j)=
     &                                               zb(index)												
                                             enddo


                                       
c>>   after reaction update mol_status
                                       
                                             do k=1,mol_tot_num
                                                mol_status_new(k)=0
                                                do k2=1,npart
                                                   if(site_status_new
     &                                                  (k2,k)
     &                                                  .eq.1)then
                                                      mol_status_new(k)
     &                                                     =1
                                                   endif
                                                enddo
                                             enddo
                                             

                                          endif
                                       endif
				
                                       
                                    endif	   
                                    
                                 endif	 
                                 
                                 
                              endif
                              
                              
                           endif
                           
                        endif
                     enddo
                  enddo
                  
               enddo
            enddo
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  2)  dimer dissociate into two subunits


            do i=1,mol_tot_num
               do j=1,npart
                  if(site_status_new(j,i).eq.1)then
                     selected_complex=0
                     selected_mol1=0
                     selected_mol2=0
                     selected_site1=0
                     selected_site2=0
                     do k=1,complex_num_new
                        if((complex_site_idx_new(k,1).eq.j)
     &                       .AND.(complex_mol_idx_new(k,1).eq.i))then
                           selected_complex=k
                           selected_mol1=complex_mol_idx_new(k,1)
                           selected_mol2=complex_mol_idx_new(k,2)
                           selected_site1=complex_site_idx_new(k,1)
                           selected_site2=complex_site_idx_new(k,2)
                        endif
                     enddo
                     if(selected_complex.ne.0)then
                        index=0
                        if((selected_site1.eq.3).AND.
     &                       (selected_site2.eq.4))then
                           index=1
                        elseif((selected_site1.eq.4).AND.
     &                          (selected_site2.eq.3))then
                           index=1
                        elseif((selected_site1.eq.5).AND.
     &                          (selected_site2.eq.5))then
                           index=2
                        endif
                        Diss_Rate_Value=Diss_Rate_max
     &                       (mol_type_index(selected_mol1),
     &                       mol_type_index(selected_mol2),
     &                       index)
                        Prob_Diss=Diss_Rate_Value*time_step
                        if(rand3(r3).lt.Prob_Diss)then
                           site_status_new(selected_site1,
     &                          selected_mol1)=0
                           site_status_new(selected_site2,
     &                          selected_mol2)=0
                           do k=selected_complex+1,complex_num_new
                              complex_mol_idx_new(k-1,1)=
     &                             complex_mol_idx_new(k,1)
                              complex_mol_idx_new(k-1,2)=
     &                             complex_mol_idx_new(k,2)
                              complex_site_idx_new(k-1,1)=
     &                             complex_site_idx_new(k,1)
                              complex_site_idx_new(k-1,2)=
     &                             complex_site_idx_new(k,2)
                           enddo
                           complex_num_new=
     &                          complex_num_new-1
                        endif
                     endif
                  endif
               enddo
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc         simulation update
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 1977       continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc   update the coordinates for all molecules

            do i=1,mol_tot_num
               do j=1,npart
                  mol_x(j,i)=mol_x_new(j,i)
                  mol_y(j,i)=mol_y_new(j,i)
                  mol_z(j,i)=mol_z_new(j,i)
                  mol_x0(j,i)=mol_x0_new(j,i)
                  mol_y0(j,i)=mol_y0_new(j,i)
                  mol_z0(j,i)=mol_z0_new(j,i)				  
                  site_status(j,i)=site_status_new(j,i)
               enddo
               mol_status(i)=0
               do j=1,npart
                  if(site_status(j,i).eq.1)then
                     mol_status(i)=1
                  endif
               enddo
            enddo
            complex_num=complex_num_new
            do i=1,complex_num
               do j=1,2
                  complex_mol_idx(i,j)=complex_mol_idx_new(i,j)
                  complex_site_idx(i,j)=complex_site_idx_new(i,j)
               enddo
            enddo

cccccccccccccc

            oligomer_num=0
            do i=1,mol_tot_num
               oligomer_subunit_num(i)=0
               do j=1,mol_tot_num
                  oligomer_subunit_index(i,j)=0
               enddo
            enddo
            
            do i=1,mol_tot_num    
               oligomer_state(i)=0        
               do j=1,mol_tot_num
                  contact_matrix(i,j)=0
                  oligomer_map(i,j)=0        
               enddo
            enddo
            do i=1,complex_num
               contact_matrix(complex_mol_idx(i,2),
     &              complex_mol_idx(i,1))=1
            enddo
            
            do i=1,mol_tot_num   
               count_flag=0
               do j=1,complex_num
                  if((complex_mol_idx(j,1).eq.i).OR.
     &                 (complex_mol_idx(j,2).eq.i))then
                     count_flag=1
                  endif
               enddo
               if(count_flag.eq.1)then
                  oligomer_state(i)=1
                  do j=1,mol_tot_num
                     color_flag(j)=0
                  enddo
                  color_flag(i)=1
                  do j=1,complex_num
                     do k=1,complex_num
                        if((color_flag(complex_mol_idx(k,1)).eq.1).AND.
     &                       (color_flag(complex_mol_idx(k,2)).eq.0)
     &                       )then
                           color_flag(complex_mol_idx(k,2))=1
                        endif
                        if((color_flag(complex_mol_idx(k,2)).eq.1).AND.
     &                       (color_flag(complex_mol_idx(k,1)).eq.0)
     &                       )then
                           color_flag(complex_mol_idx(k,1))=1
                        endif
                     enddo
                  enddo
                  do j=1,mol_tot_num
                     oligomer_map(i,j)=color_flag(j)
                  enddo
               endif
            enddo
            
            do i=1,mol_tot_num
               if(oligomer_state(i).eq.1)then
                  count_flag=0
                  do j=1,oligomer_num
                     do k=1,oligomer_subunit_num(j)
                        if(oligomer_subunit_index(j,k).eq.i)then
                           count_flag=1
                        endif
                     enddo
                  enddo
                  if(count_flag.eq.0)then
                     oligomer_num=oligomer_num+1
                     do j=1,mol_tot_num
                        if(oligomer_map(i,j).eq.1)then
                           oligomer_subunit_num(oligomer_num)=
     &                          oligomer_subunit_num(oligomer_num)+1
                           oligomer_subunit_index(oligomer_num,
     &                          oligomer_subunit_num(oligomer_num))=j
                        endif
                     enddo
                  endif
               endif
            enddo
            
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c  output data along the trajectory
ccccccccccccccccccccccccccccccccccccccccccccccccccc    
c            print*,n_t,mc_time_step,complex_num
c            do i=1,complex_num
c               print*,'interface',i,complex_mol_idx(i,1),
c     &              complex_site_idx(i,1),complex_mol_idx(i,2),
c     &              complex_site_idx(i,2)
c            enddo

            if((output_flag_rec.eq.1).AND.
     &           (MOD(mc_time_step,rec_output_freq).eq.0))then    
               open (unit=10,file=
     &              'CageAssAlignFlex_rec_VL100.dat',
     &              status='unknown',access='append')
               write(10,2200) 'index',n_t,mc_time_step,complex_num,
     &              oligomer_num
               do i=1,oligomer_num
                  write(10,2500) i,'oligomer size: ',
     &                 oligomer_subunit_num(i)
               enddo
               close(10)
            endif

            if((output_flag_str.eq.1).AND.
     &           (MOD(mc_time_step,trj_output_freq).eq.0))then    
               open (unit=10,file=
     &              'CageAssAlignFlex_traj_VL100.pdb'
     &              ,status='unknown',access='append')
               write(10,2200) 'index',n_t,mc_time_step,complex_num
               do i=1,mol_tot_num
                  if(mol_type_index(i).eq.1)then
                     write(10,2100)'ATOM  ',1,' CA ','ALA', 
     &                    'A',i,mol_x(1,i),mol_y(1,i),
     &                    mol_z(1,i)
                     write(10,2100)'ATOM  ',2,' CA ','GLY', 
     &                    'A',i,mol_x(2,i),mol_y(2,i),
     &                    mol_z(2,i)
                     write(10,2100)'ATOM  ',3,' CA ','LEU', 
     &                    'A',i,mol_x(3,i),mol_y(3,i),
     &                    mol_z(3,i)
                     write(10,2100)'ATOM  ',4,' CA ','LEU', 
     &                    'A',i,mol_x(4,i),mol_y(4,i),
     &                    mol_z(4,i)
                     write(10,2100)'ATOM  ',5,' CA ','ILE', 
     &                    'A',i,mol_x(5,i),mol_y(5,i),
     &                    mol_z(5,i)
                     write(10,2100)'ATOM  ',6,' CA ','PRO', 
     &                    'A',i,mol_x(6,i),mol_y(6,i),
     &                    mol_z(6,i)
                     write(10,2100)'ATOM  ',7,' CA ','PRO', 
     &                    'A',i,mol_x(7,i),mol_y(7,i),
     &                    mol_z(7,i)
                     write(10,2100)'ATOM  ',8,' CA ','VAL', 
     &                    'A',i,mol_x(8,i),mol_y(8,i),
     &                    mol_z(8,i)	 
                     write(10,3000) 'TER'
                  endif
               enddo
               write(10,3000) 'END'      
               close(10)

            endif

 2100       format(A6,I5,1x,A4,1x,A3,1x,A1,I4,4x,3F8.3)
 2200       format(A5,I5,I10,I5,1x,I5)
 2300       format(I3,4I5)
 2400       format(I3,I3,I5)
 2500       format(I3,1x,A15,I5)
 3000       format(A3)     

ccccccccccccccccccccccccccccccccccc
cc  end current simulation step
ccccccccccccccccccccccccccccccccccc

         enddo




cccccccccccccccccccccccccccc
cc end current trajectory
cccccccccccccccccccccccccccc

      enddo

cccccccccccccccccccc

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc  
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                                  
ccccccccccccccccccccccccccccccccccccccccccccccccc            
   
ccccccccccccccccccccccccccccccccccccccccccccccccc            
            
      subroutine gettheta(point_x,point_y,point_z,test_theta)

      implicit none
      real*8 point_x(3),point_y(3),point_z(3),test_theta
ccc
      real*8 lx(2),ly(2),lz(2),lr(2)
      real*8 conv,doth1,doth2,conv2
      integer in

ccccccccccccccccccc
c>>  bond lengt
ccccccccccccccccccc

      do in=1,2
         lx(in)=0
         ly(in)=0
         lz(in)=0
         lr(in)=0
      enddo

      lx(1)=point_x(2)-point_x(1)
      ly(1)=point_y(2)-point_y(1)
      lz(1)=point_z(2)-point_z(1)
      lr(1)=sqrt(lx(1)**2+ly(1)**2+lz(1)**2)
      
      lx(2)=point_x(3)-point_x(2)
      ly(2)=point_y(3)-point_y(2)
      lz(2)=point_z(3)-point_z(2)
      lr(2)=sqrt(lx(2)**2+ly(2)**2+lz(2)**2)
 
cccccccccccccccccccc
c>>  theta value
cccccccccccccccccccc
      
      test_theta=0

      conv=180/3.14159 
      do in=1,1
         doth1=-1*(lx(in+1)*lx(in)+ly(in+1)*ly(in)+lz(in+1)*lz(in))
         doth2=doth1/(lr(in+1)*lr(in))
         if(doth2.gt.1)then
            doth2=1
         endif
         if(doth2.lt.-1)then
            doth2=-1
         endif
         test_theta=acos(doth2)*conv
      enddo

cccccccccccccccc

      return
      end

cccccccccccccccccccccccccccccccccccc
                        
