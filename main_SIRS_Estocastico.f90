
module sirs_estocastico
   use geraRede
   use mod_rndgen
   use mod_tools
   implicit none
   
   integer, parameter :: dp = kind(0.0d0)
!#######################################################################      
    real(dp), allocatable :: Pn_QS_global(:)
	real(dp):: rho_medioQS_global
	real(dp) :: rho2_medioQS_global
	real(dp) :: dev_rho_medioQS_global 
	real(dp) :: Xi_global
	real(dp) :: S_Schanon_global
	real(dp) :: t_LS
	real(dp) :: Y4
	real(dp), allocatable :: v1(:)
	integer :: ind_v1
	real(dp) :: qtdade
!#######################################################################      
!   real(dp), allocatable :: I_i(:)
!   real(dp), allocatable :: R_i(:)
!#######################################################################
!  Tem tamanho this%list e guarda os estados 0, 1, 2 para cada sitio         
!   integer(kind=1), allocatable :: sigma(:)
!#######################################################################
!  nInf eh o numero de infectados, enquanto nRec eh o numero de
!  recuperados. Quando um sítio i1 eh infectado, nInf = nInf + 1,
!  sigma(i1) = 1 e o novo sitio infectado entra na lista inf_List:
!  inf_List(nInf) = i1.
!  Quando um sitio na posicao i1 sorteado da lista inf_List, infectado,
!  se recupera, nRec = nRec + 1 e rec_List(nRec) = inf_List(i1),
!  sigma(rec_List(nRec)) = 2, 
!  inf_List(i1) = inf_List(nInf), nInf = nInf - 1 (nesta ordem). 
!  Quando um sitio recuperado, escolhido ao acaso da rec_List, 
!  rec_List(i1), torna-se suscetivel, rec_List(i1) = rec_List(nRec),
!  sigma(rec_List(i1)) = 1 e nRec = nRec - 1.
!  De forma analoga com nRec e rec_List.
   integer :: nInf
   integer :: nRec
   integer :: nSus   
   integer :: stubs_inf, stubs_infTotal
   integer, allocatable :: inf_List(:)
   integer, allocatable :: rec_List(:)
   integer, allocatable :: sus_List(:)
   real(dp), allocatable :: instante_da_foto(:)    
   logical :: bandeira
   real(dp) :: valor_subcritico
   real(dp) :: t_teste
!#######################################################################
!  n_ens eh o numero de copias que guardamos para backup da dinamica
!  quando o sistema alcanca o estado absorvente.
!  Abaixo do ponto critico se fara muito necessario esse backup,
!  uma vez que o sistema visita bastante o estado absorvente.
!  Acima do ponto critico isso raramente acontece, a nao ser por alguma
!  flutuacao, embora para tempos suficientemente grandes, a dinamica
!  eventualmente visita o estado absorvente.
!  Por sua vez, o residuo eh o numero maximo de sitios infectados
!  e recuperados que guardaremos em cada uma dessas copias de backup.
!  ens_Inf(:,:) eh uma matriz na qual cada coluna representa
!  uma configuracao e cada linha representa sitios no estado infectado
!  naquela configuração.
!  Analogamente para a matriz ens_Rec(:,:).
!  Note, n_ens eh o numero maximo!
!  Em ens_nInf(:), um vetor, guardamos o numero de sitios infectados
!  naquela configuracao. A ens_nInf(i1), membro i1 da colecao,
!  corresponde uma confugracao ens_Inf(1:ens_nInf(i1), i1)
!  com ens_nInf(i1) sitios infectados.
!  O analogo ocorre com ens_Rec(:,:).
  
   integer :: n_ens
   integer :: residuo
   integer, allocatable :: ens_Inf(:,:)
   integer, allocatable :: ens_Rec(:,:)
   integer, allocatable :: ens_nInf(:)
   integer, allocatable :: ens_nRec(:)
!#######################################################################         
   real(dp) :: Im, Rm, Sm
   real(dp) :: t, dt, t_relax, t_max
!#######################################################################      
   integer :: nInfMax
   !character(len=500), parameter :: local = trim(adjustl("/home/jota/SIRS_Estocastico/Rst_Rede_Real/"))
   integer(kind=8) :: n_ite
   
!#######################################################################
   contains

!#######################################################################    
   subroutine aloca_listas_dinamicas(this)
      class(grafo), intent(in) :: this
      integer :: j1
      if(.not. allocated(lista_de_clusters))then
         write(*,*) ""
         write(*,*) "Eh preciso rodar a rotina: "
         write(*,*) ""
         write(*,*) "sub_classifica_clusters,"
         write(*,*) ""
         write(*,*) "do modulo mod_tools_redes.f90"
         write(*,*) ""
         stop
      endif
!#######################################################################
      if(allocated(Pn_QS_global)) deallocate(Pn_QS_global)
         allocate(Pn_QS_global(comp_gigante))    
!#######################################################################
      if(allocated(sigma)) deallocate(sigma)
         allocate(sigma(this%nodes))
!#######################################################################
      if(allocated(inf_List)) deallocate(inf_List)
         allocate(inf_List(comp_gigante))
!#######################################################################
      if(allocated(rec_List)) deallocate(rec_List)
         allocate(rec_List(comp_gigante))
!#######################################################################
!      if(allocated(sus_List)) deallocate(sus_List)
!         allocate(sus_List(comp_gigante))  
!#######################################################################
!      if(allocated(instante_da_foto)) deallocate(instante_da_foto)
!         allocate(instante_da_foto(comp_gigante))           
!#######################################################################
         n_ens = 100
         residuo = int( qtdade * comp_gigante)
!#######################################################################
   if(allocated(ens_Inf)) deallocate(ens_Inf)
      allocate(ens_Inf(residuo, n_ens))
!#######################################################################
   if(allocated(ens_nInf)) deallocate(ens_nInf)
      allocate(ens_nInf(n_ens))
!#######################################################################
   if(allocated(ens_Rec)) deallocate(ens_Rec)
      allocate(ens_Rec(residuo, n_ens))
!#######################################################################
   if(allocated(ens_nRec)) deallocate(ens_nRec)
      allocate(ens_nRec(n_ens))
!#######################################################################                             

      stubs_infTotal = 0
      do j1 = 1, this%nodes
         if(lista_de_clusters(j1) /= i_comp_gigante)cycle
         stubs_infTotal = stubs_infTotal + this%deg(j1)
      enddo
      
      if(allocated(v1)) deallocate(v1)
         allocate(v1(this%nodes))
                  
   end subroutine
   
!#######################################################################   
   
   subroutine condicao_inicial(this, alp, lamb, mu, t_0, tMax, tRelax)
!     integer, parameter :: dp = kind(0.0d0)
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp, lamb, mu
      real(dp), intent(in) :: t_0, tMax, tRelax
      integer :: j1, l10
!#######################################################################

      nInf = comp_gigante
      nInfMax = comp_gigante
      nRec = 0
      nSus = 0
      
      if(.not.allocated(lista_de_clusters)) stop "Rode a subrotina sub_classifica_clusters()"
!#######################################################################
      l10 = 0
      do j1 = 1, this%nodes
         if(lista_de_clusters(j1) /= i_comp_gigante)then
            sigma(j1) = 3
            cycle
         endif
         l10 = l10 + 1
         inf_List(l10) = j1
         sigma(j1) = 1 
      enddo
      stubs_inf = stubs_infTotal
      !write(*,*) "stubs_infTotal = ", stubs_infTotal
!#######################################################################     
      !rec_List = 0
      Pn_QS_global = 0.0_dp
 !######################################################################     
      v1 = 0.0_dp
               
      !ens_Inf = 0
      !ens_nInf = 0
      !ens_Rec = 0
      !ens_nRec = 0
      
      t = t_0
      t_max = tMax
      t_relax = tRelax 
 !######################################################################            
   end subroutine
!#######################################################################   
   subroutine sirs_estoc(this, alp, lamb, mu, T_vs)
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp, lamb, mu
      logical, intent(in) :: T_vs
      real(dp) :: Tax_inf, Tax_rec, Tax_wan
      real(dp) :: Tax_total
      real(dp) :: prob_inf, prob_rec, prob_wan
      logical :: recup, recai, infec
      real(dp) :: prob
      type(rndgen) :: gen
      integer :: seed_est
      integer :: ultima_foto

      integer, parameter :: nfotos = 100
      integer :: j1, j2, j3,j4, k1, k2, k3, k4
      real(dp) :: probInf, prob_dist_t
      real(dp), parameter :: p_deco = 0.02_dp
      real(dp) :: prob1, p_act, p_fill
      integer :: n_vezes
      integer :: sumPQS      
      character(len=100) :: dados
      character(len=10) :: alp_char, lamb_char
      character(len=50) :: rede_char
      
      seed_est = 999999991
      call gen%init(seed_est)
!#######################################################################
      if( .not. allocated(ens_Inf)) stop "Provavelmente a rotina aloca_listas_dinamicas &
      nao foi executada. Faca isso."
!#######################################################################      
      ultima_foto = 0

!#######################################################################       
ld:   do while(t <= t_max)
!#######################################################################         
         Tax_inf = 1.0_dp * stubs_inf * lamb
         Tax_rec = 1.0_dp * nInf * mu
         Tax_wan = 1.0_dp * nRec * alp     
!#######################################################################             
         Tax_total = Tax_inf + Tax_rec + Tax_wan
!#######################################################################
         prob_rec = 1.0_dp * Tax_rec/Tax_total
         prob_inf = 1.0_dp * Tax_inf/Tax_total
         prob_wan = 1.0_dp * Tax_wan/Tax_total
!#######################################################################        
         prob = gen%rnd()
         dt = -1.0_dp * log(max(1e-12,prob))/Tax_total
         t = t + dt
!#######################################################################
         if(t >= t_relax) Pn_QS_global(nInf) = Pn_QS_global(nInf) + dt
!#######################################################################
         if(t > 5_dp)then
!#######################################################################
            if(ultima_foto < nfotos)then
               prob1 = gen%rnd()
               p_fill = 1_dp * dt
!#######################################################################               
               if(prob1 <= p_fill)then
                  
                     ultima_foto = ultima_foto + 1
                     ens_nInf(ultima_foto) = nInf
                     ens_nRec(ultima_foto) = nRec
!#######################################################################                 
                     do j1 = 1, nInf
                        ens_Inf(j1, ultima_foto) = inf_List(j1)                     
                     enddo
                     !write(*,*) nInf, " no tempo", t
!#######################################################################                  
                     do j1 = 1, nRec                     
                        ens_Rec(j1, ultima_foto) = rec_List(j1)
                     enddo                
!#######################################################################
               endif
!#######################################################################
            else
               prob1 = gen%rnd()
               p_act = 1.0_dp * p_deco * dt
!#######################################################################               					
               if(prob1 <= p_act)then
                     j1 = gen%int(1, ultima_foto)
                     ens_nInf(j1) = nInf
                     ens_nRec(j1) = nRec       
                     do j2 = 1, nInf
                        ens_Inf(j2, j1) = inf_List(j2)
                     enddo
                     do j2 = 1, nRec
                        ens_Rec(j2, j1) = rec_List(j2)                   
                     enddo
               endif
!#######################################################################
            endif						
!#######################################################################
         endif
!#######################################################################
 ! Aqui acontecem as transicoes
!#######################################################################
       prob = gen%rnd()
      
       if(prob <= prob_rec)then
            if(t < 20.0_dp)then
!#######################################################################
! Se t < 20.0, usa-se a condicao de contorno refletora
!#######################################################################              
               if(nInf == 1)then
                  
                  sigma(inf_List(nInf)) = 2
                  nRec = nRec + 1
                  rec_List(nRec) = inf_List(nInf)
                  stubs_inf = 0
                  prob = gen%rnd()

                  if(prob <= 1.0_dp * nRec/comp_gigante)then                     
                     j3 = gen%int(1, nRec)
                     j1 = rec_List(j3)
                     rec_List(j3) = rec_List(nRec)
                     nRec = nRec - 1
                  else
lant:                do
                       j1 = gen%int(1, this%nodes)
                       if(lista_de_clusters(j1) /= i_comp_gigante) cycle lant
                       if(sigma(j1) /= 0)cycle lant
                       exit lant
                     enddo lant
                  endif

                  inf_List(nInf) = j1
                  sigma(j1) = 1
                  stubs_inf = this%deg(j1)
                  cycle ld			
               endif
            endif         
               
            j3 = gen%int(1, nInf)
            j1 = inf_List(j3)
            
            if(t >= t_relax) v1(j1) = v1(j1) + 1.0_dp/mu
            
            sigma(j1) = 2
            nRec = nRec + 1
            rec_List(nRec) = j1

            inf_List(j3) = inf_List(nInf)
            nInf = nInf - 1
            stubs_inf = stubs_inf - this%deg(j1)
!#######################################################################################################################
!		Se o estado absorvente eh alcancado, escolhemos uma configuracao ja visitada

			if(nInf == 0)then
                                
                                do j3 = 1, nRec
                                   sigma(rec_List(j3)) = 0
                                enddo
                                
				k1 = gen%int(1,ultima_foto)

				nInf = ens_nInf(k1)
				nRec = ens_nRec(k1)
				
				do j3 = 1, nInf
                                     
				   inf_List(j3) = ens_Inf(j3, k1)

				   j1 = inf_List(j3)
				   sigma(j1) = 1	
				   stubs_inf = stubs_inf + this%deg(j1)
                                enddo
!#######################################################################
!
!Aqui eh o primeiro lugar onde aparece um sitio de indice negativo!
!
!#######################################################################
				do j3 = 1, nRec
				   rec_List(j3) = ens_Rec(j3, k1)
				   j1 = rec_List(j3)
				   sigma(j1) = 2
				enddo
			endif
         elseif(prob <= (prob_rec+prob_inf))then
sel:           do
                  j3 = gen%int(1,nInf)
                  j1 = inf_List(j3)
                  prob = gen%rnd()
                  if(prob <= (1.0_dp * this%deg(j1)/this%degMax)) exit sel
               enddo sel
               k1 = gen%int(this%aux(j1), this%aux(j1) + this%deg(j1) - 1)
               k2 = this%listAdj(k1)
               if(sigma(k2) /= 0)cycle ld
               sigma(k2) = 1
               nInf = nInf + 1
               inf_List(nInf) = k2
               stubs_inf = stubs_inf + this%deg(k2)
          elseif(prob <= (prob_rec+prob_inf+prob_wan))then
               j3 = gen%int(1,nRec)
               j1 = rec_List(j3)
               sigma(j1) = 0
               rec_List(j3) = rec_List(nRec)
               nRec = nRec - 1
         endif
      enddo ld
      
      if(T_vs)then
         close(888)
      endif
	!###########################################################
	!	Serve tanto para a Pn_QS, quanto para Pn_QS_n_redundante	
		sumPQS = sum(Pn_QS_global)
	!###########################################################					
      				
	do j1 = 1, size(Pn_QS_global)			
		Pn_QS_global(j1) = 1.0_dp * Pn_QS_global(j1) / sumPQS
	enddo
							
	rho_medioQS_global = 0.0_dp
	rho2_medioQS_global = 0.0_dp
	!S_Schanon_global = 0.0_dp			
	
	if(Pn_QS_global(1) > 0.0_dp)then
		t_LS = 1.0_dp /Pn_QS_global(1)
	else
		t_LS = t_max
	endif

	
	do j1 = 1, nInfMax
		if(Pn_QS_global(j1) > 0.0_dp)then
			rho_medioQS_global = 1.0_dp * j1 * Pn_QS_global(j1) + rho_medioQS_global						! n medio
			rho2_medioQS_global =  ((1.0_dp *j1) ** 2.0_dp) * Pn_QS_global(j1) + rho2_medioQS_global
			!S_Schanon_global = S_Schanon_global - Pn_QS_global(j1) * log(Pn_QS_global(j1))
		endif
	enddo

	rho_medioQS_global = 1.0_dp * rho_medioQS_global / comp_gigante

	rho2_medioQS_global = 1.0_dp * rho2_medioQS_global/ ((1.0_dp * comp_gigante)**2.0_dp)

	dev_rho_medioQS_global = (rho2_medioQS_global - rho_medioQS_global**2.0_dp)**0.5_dp
	
	Xi_global = 1.0_dp * comp_gigante * (rho2_medioQS_global - (rho_medioQS_global**2.0_dp))/rho_medioQS_global
	
	
		
	v1 = v1/(1.0_dp * (t - t_relax))
	
	v1 = v1/(sum(v1**2.0_dp))**0.5_dp

	Y4 = sum(v1**4.0_dp)
	
   end subroutine
!#######################################################################
   
end module


module limiaresCampoMedio
   use geraRede
   use mod_rndgen
   implicit none
   
   integer, private, parameter :: dp = kind(0.0d0)
   
   real(dp), allocatable :: x(:), y(:)
   
   contains
   
   subroutine aloca_listas_e_matrizes(this)
      class(grafo) :: this
      if(allocated(x)) deallocate(x)
      allocate(x(this%nodes))

      if(allocated(y)) deallocate(y)
      allocate(y(this%nodes))
      
   end subroutine
   
   subroutine limiarQMF(this, lambda1)
      class(grafo), intent(in) :: this
      real(dp), intent(inout) :: lambda1
      real(dp), parameter :: tole = 1.0d-4
      !=================================================================
      real(dp) :: autovalor0 
      !=================================================================
      integer :: j1, j2, j3, j4
      real(dp) :: xp, yp
      integer :: ipx, ipy
      type(rndgen) :: geni
      real(dp) :: soma
      integer :: vizim
      real(dp) :: a_l, a_2l
      integer :: semente2
      real(dp) :: erro

      call aloca_listas_e_matrizes(this)
      !=================================================================
      
      !=================================================================      
      call geni%init(semente2)     
      !=================================================================
	  ! Condicao inicial				
	  !x = 1.0d0
	  !xp = 1.0d0
      !ipx = 1
      
      xp = 0.0d0
      ipx = 1
      do j1 = 1, this%nodes
         x(j1) = geni%rnd() * 1.0d-3/(1.0d0 * this%nodes)
         if( abs(x(j1)) > abs(xp) )then
            xp = x(j1)
            ipx = j1
         endif
      enddo
      x = x/xp
      xp = 1.0d0
      !=================================================================
      ! Apos normalizar o vetor x, a norma ||xp|| = ||x||_{oo}
      ! se torna = 1.0d0.
      !=================================================================
      j1 = 0
      !=================================================================  
!=======================================================================      
iter: do while(.True.)	
	     !Aqui comeca o algoritmo
         !==============================================================
         yp = 0.0d0
         ipy = 1
         
         !==============================================================
         ! Y = A X
         !==============================================================
         do j2 = 1, this%nodes
            !===========================================================
            soma = 0.0_dp
            do j3 = this%aux(j2), this%aux(j2) + this%deg(j2) - 1
               vizim = this%listAdj(j3)
               soma = soma + x(vizim)
            enddo
            y(j2) = soma
            !===========================================================
            if( abs(y(j2)) > abs(yp) )then
               yp = y(j2)
               !========================================================
               ! Esse ipy sera o proximo ipx.
               !========================================================
               ipy = j2
            endif
            !===========================================================
         enddo
         !write(*,*) "y"
         !write(*,*) y
         !stop
         !==============================================================
         ! A componente y(ipx) representa o ganho multiplicativo
         ! obtido com o produto matricial.
         !==============================================================
         autovalor0 = y(ipx)
		 !==============================================================
		 ! Se yp = 0, para tudo!
	     if(yp == 0.0d0)then
	        xp = 0.0d0
	        ipx = 1
            do j2 = 1, this%nodes   
			   x(j2) = 2.0d0 * geni%rnd()
			   if( abs(x(j2)) > abs(xp) )then
			      xp = x(j2)
			      ipx = j2
			   endif
			enddo
			x = x/xp
			xp = 1.0d0
			!===========================================================
			! Como o vetor x foi normalizado via || ||_{oo},
			! xp = 1.0d0
			!===========================================================

			j1 = 0
			write(*,*) 'Autovalor nulo encontrado. Retomando outro vetor x...'
			cycle iter
         endif
         !==============================================================
		 erro = abs(maxval(x - y/yp))		 
		 !==============================================================
         x = y/yp
         !==============================================================
         !write(*,*) "==================================================="
         !write(*,*) "Erro = ", erro, " tole = ", tole 
         !write(*,*) "==================================================="
         !call sleep(5)
          
         
         if(erro < tole)then
            write(*,*) "LEV Jacobiana = ", autovalor0
            lambda1 = 1.0_dp/autovalor0
            write(*,*) "lambda_C = ", 1.0_dp/autovalor0            
            exit iter
         endif
         
         !==============================================================
		 ! O ipy achado apos o produto matricial eh o indice ipx da
		 ! primeira maior componente do vetor x.
		 ! xp neste caso agora = 1.0d0, pois y ja foi normalizado
		 ! na linha anterior ao comando 'if(erro < tole)then'
		 ipx = ipy
		 xp = 1.0d0
		 !==============================================================
         j1 = j1 + 1
         !==============================================================
      enddo iter   
      
      deallocate(x)
      deallocate(y)     
   end subroutine
   
end module


program main
   use geraRede
   use sirs_estocastico
   use mod_rndgen
   use mod_tools
   use limiaresCampoMedio
   implicit none
!#######################################################################   
   type(grafo_PL_UCM) :: rede
!#######################################################################   
   real(dp) :: dlamb, dlamb2
   real(dp) :: lamb0
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
!#######################################################################
   !type(rndgen) :: gen1
   integer :: seed(10)
   integer :: seed1
   type(rndgen) :: ger_inic
!#######################################################################
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
!#######################################################################   
   character(len=500) :: t_vs_Im
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo
   character(len=5) :: tam_char
   character(len=5) :: gama_char
   character(len=5) :: indice
!#######################################################################   
   integer :: i1, i2, i3
   integer :: ind_amostra
   logical :: T_vs
!#######################################################################
   integer :: sumDeg2
   
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   real(dp) :: tTotal, tRex
   real(dp) :: lambdaf
   real(dp) :: qm, q2m
   integer :: nargus
   character(len=20) :: buffer
   integer :: divisor
   character(len=100) :: cwd, resultados, tipoCorte
   character(len=400) :: local
   real(dp) :: rho_aux, Xi_aux, Xi_max
   real(dp) :: lamb1, lamb2, lambC
   logical :: existe
   integer :: status_io
   character(len=10) :: interpola
   logical :: interpola_bool
   logical :: TeveLeitura
   integer :: ind_interp
   integer :: int_soCalculaIPR
   logical :: soCalculaIPR
   logical :: taAberta
   real(dp) :: sumIi2   
   intger :: gasta_aleatorio
   !=================================================
   call getcwd(cwd)
   local = trim(adjustl(cwd))//'/'

   !tipoCorte ='_Rigido'
   !tipoCorte ='_2sqrtN'
   tipoCorte ='_sqrtN'
   
   resultados = 'Rst_Estocastico_Corte'//trim(adjustl(tipoCorte))
   
   resultados = trim(adjustl(resultados))

   call system('mkdir -p '//resultados)
    
   local = trim(adjustl(local))//trim(adjustl(resultados))//"/"
   
   
!#######################################################################   
   seed1=947361823
!#######################################################################   

   nargus = iargc()

   if(nargus == 14)then
      !#############################
      ! Amostra
      !#############################
      call getarg(1, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) ind_amostra
      !#############################
      write(indice,'(I0)') ind_amostra
      !#############################
      !	Tamanho da rede
      !#############################
      call getarg(2, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) tam_rede

      !#############################
      ! Tamanho	da rede
      !#############################
      call getarg(3, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) grau_min

      !#############################
      ! Expoente Gama
      !#############################
      call getarg(4, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) gama_exp

      !#############################
      ! Lambda0
      !#############################
      call getarg(5, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) lamb0

      !#############################
      ! Divisor que fornece dlambda
      !#############################
      call getarg(6, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) divisor
      
      write(*,*) "O valor do divisor de 0.0125 eh: ", divisor

      !#############################
      ! Lambdaf
      !#############################
      call getarg(7, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) lambdaf

      !#############################
      ! Alfa
      !#############################
      call getarg(8, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) alp

      !#############################
      ! Tempo relaxacao
      !#############################
      call getarg(9, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) tRex

      !#############################
      ! Alfa
      !#############################
      call getarg(10, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) tTotal
      
      call getarg(11, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) interpola

      if( trim(adjustl(interpola)) == 'false' )then
         interpola_bool = .False.
         write(*,*) 'main_SIRS_Estocastico nao interpolara'
      elseif( trim(adjustl(interpola)) == 'true' )then
         interpola_bool = .True.
         write(*,*) 'main_SIRS_Estocastico interpolara'
      else
         stop 'Valor de interpola nao pode ser estabelecido'
      endif
      
      call getarg(12, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) ind_interp
      
      call getarg(13, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) int_soCalculaIPR
      
      if(int_soCalculaIPR == 0)then
         soCalculaIPR = .False.
      elseif(int_soCalculaIPR == 1)then
         soCalculaIPR = .True.
      else
         stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
      endif

      call getarg(14, buffer)
      buffer = trim(adjustl(buffer))
      read(buffer,*) qtdade

   else
      stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
   endif
!#######################################################################

   dlamb = 0.0125_dp/(1.0_dp * divisor)

   write(*,*) "######################################"
   write(*,*) "t Total = ", int(tTotal)
   write(*,*) "t relaxamento = ", int(tRex)
!#######################################################################
   write(*,*) "########Parametros iniciais e simulacao#########"
   write(*,*) "dlamb = ", dlamb
!#######################################################################   
   if(tam_rede == 1000)then
      tam_char = '1k'
   elseif(tam_rede == 3000)then
      tam_char = '3k'
   elseif(tam_rede == 10000)then
      tam_char = '10k'
   elseif(tam_rede == 30000)then
      tam_char = '30k'
   elseif(tam_rede == 100000)then
      tam_char = '100k'
   elseif(tam_rede == 300000)then
      tam_char = '300k'
   elseif(tam_rede == 1000000)then
      tam_char = '1M'
   elseif(tam_rede == 3000000)then
      tam_char = '3M'
   elseif(tam_rede == 10000000)then
      tam_char = '10M'
   elseif(tam_rede == 30000000)then
      tam_char = '30M'
   elseif(tam_rede == 100000000)then
      tam_char = '100M'      
   else
      stop 'Escolha um tamanho de rede dentro do catalogo'
   endif
!#######################################################################
   if( trim(adjustl(resultados)) == 'Rst_Estocastico_Corte_Rigido' )then 
      call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
   elseif(trim(adjustl(resultados)) == 'Rst_Estocastico_Corte_sqrtN' )then
      grau_max = (1.0_dp * tam_rede)**(0.5_dp)
   elseif(trim(adjustl(resultados)) == 'Rst_Estocastico_Corte_2sqrtN' )then
      grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
   endif   
!#######################################################################   
   write(*,*) "##############Grau Maximo Estimado####################"    
   write(*,*) "O grau maximo estimado para a rede eh: ", grau_max
!#######################################################################
   write(gama_char,'(f5.2)') gama_exp
!#######################################################################
   call ger_inic%init(seed1)
   i2 = 1
   do i1 = 1, 1000
      if(mod(i1,100) > 0)then
         gasta_aleatorio = ger_inic%int(100000000,999999999)
         cycle
      endif
      seed(i2)  = ger_inic%int(100000000,999999999)
      write(*,*) i1, seed(i2)
      i2 = i2+1      
   enddo

   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
   
   local = trim(adjustl(trim(adjustl(local))//'gam_'//trim(adjustl(gama_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )   

   local = trim(adjustl(trim(adjustl(local))//'ams_'//trim(adjustl(indice))//'/'))
   
   call system('mkdir -p '//trim(adjustl(local)) ) 

!=======================================================================
!				Inicia grafo
!=======================================================================

   call criaRedeEClassificaClusters(rede, tam_rede, grau_min, grau_max, gama_exp, seed(ind_amostra))
!=======================================================================
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!=======================================================================


   write(*,*) "#########Concluiu Ligacao da Rede##########"
   write(*,*) "Tamanho da rede eh ", rede%nodes
   write(*,*) ""
   write(*,*)"O grau minimo da rede eh: ", rede%degMin
   write(*,*) ""
   write(*,*)"O grau maximo da rede eh: ", rede%degMax
   write(*,*) ""
   write(*,*)"O grau medio da rede eh: ", rede%degMean
   write(*,*) ""
   write(*,*) "O expoente da distribuicao PL eh ", gama_exp
   write(*,*) "" 

   write(alp_char2, '(f9.3)') alp   

   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//trim(adjustl(tipoCorte))//'/'))
      

   call system('mkdir -p '//trim(adjustl(local)) ) 
   
   !--------------------------------------------------------------------
   call kNN_e_clustering(rede)
   !--------------------------------------------------------------------
   call calcula_Pk_E_Lambda0(rede)
   !--------------------------------------------------------------------   
!#######################################################################   
   call aloca_listas_dinamicas(rede)
!#######################################################################
!#######################################################################
   write(*,*) "######Parametros Dinamicos da Epidemia######"   
   write(*,*) "O valor de alp eh ", alp
   write(*,*) "O valor de mu eh ", mu
!#######################################################################
  

   if(ind_interp > 0)then
      arquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_XiEst.dat'))
      inquire(file=trim(adjustl(arquivo)), exist=existe)      
      if(existe)then
         write(*,*) 'Aquivo existe e esta sendo lido'
         open(333, file=trim(adjustl(arquivo)), status='old')
         read(333, *, iostat=status_io) lamb1, Xi_aux
         Xi_max = Xi_aux !Provisorio
         if(status_io == 0)then
            lamb0 = lamb1
            !write(*,*) 'Lambda0 da interpolacao foi lido'
         else
            stop 'lambda0 da interpolacao nao pode ser lido. Interpolacao abortada.'
         endif
         read(333, *, iostat=status_io) lamb2, Xi_aux
         
         if(status_io == 0)then
            dlamb2 = lamb2 - lamb1
            !write(*,*) 'Lambda2 foi lido = ', lamb2
         else
            stop 'Lambda2 nao foi lido. Interpolacao abortada.'
         endif
         rewind(333)
                  
         do
            read(333, *, iostat=status_io) lamb1, Xi_aux
            if(status_io == 0)then
               if( Xi_aux > Xi_max)then
                  lamb0 = lamb1
                  Xi_max = Xi_aux
               endif
            else
               exit
            endif
         enddo

         lamb0 = lamb0 - 3.0d0 * dlamb
         nlamb = 12
         lambdaf = lamb0 + 1.0_dp * nlamb * dlamb

         if(ind_interp == 1)then
            lamb0 = lamb0 + 1.0d0 * dlamb/4.0d0
         elseif(ind_interp == 2)then
            lamb0 = lamb0 + 1.0d0 * dlamb/3.0d0
         elseif(ind_interp == 3)then
            lamb0 = lamb0 + 1.0d0 * dlamb/2.0d0
         elseif(ind_interp == 4)then
            lamb0 = lamb0 + 2.0d0 * dlamb/3.0d0
         elseif(ind_interp == 5)then
            lamb0 = lamb0 + 3.0d0 * dlamb/4.0d0
         elseif(ind_interp == 6)then
            lamb0 = lamb0
         else
            stop "Entre com valores de 1 a 6. Abortando programa..."
         endif
         write(*,*) "Ind_interp = ", ind_interp
         close(333)

      endif
   elseif(ind_interp == 0)then
      arquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rhoEst.dat'))
      inquire(file=trim(adjustl(arquivo)), exist=existe)
      if(existe)then
         TeveLeitura = .False.
         open(333, file = trim(adjustl(arquivo)), status = 'old')
llbd0:   do
            read(333, *, iostat = status_io) lamb0, rho_aux
            if(status_io /= 0)then
               exit llbd0
            else
               TeveLeitura = .True.   
            endif
         enddo llbd0
         close(333)
         lamb0 = lamb0 + dlamb
      endif
   endif
   
   nlamb = int(((lambdaf - lamb0)/dlamb))
   !--------------------------------------------------
   write(*,*) "n_lamb = ", nlamb
   write(*,*) "lambda0 = ", lamb0
   write(*,*) "lambdaf = ", lambdaf
   !--------------------------------------------------
   lamb = lamb0
   !--------------------------------------------------
   write(*,*) 'O valor de lamb0 = ', lamb
      
   if(soCalculaIPR)then
      inquire(file=trim(adjustl(local))//trim(adjustl('lbd_vs_XiEst.dat')), exist=existe)
      
      if(existe)then
         open(334, file=trim(adjustl(local))//trim(adjustl('lbd_vs_XiEst.dat')), status='old')
      else
         stop "Nao eh possivel achar o IPR no limiar, pois nao existe arquivo lbd_vs_XiEst.dat"
      endif
      
      Xi_max = 0.0_dp
      
      do
         read(334, *, iostat=status_io) lamb, Xi_global
         
         if(status_io == 0)then
            if(Xi_global > Xi_max)then
               Xi_max = Xi_global
               lambC = lamb
            endif
         else
            exit
         endif         
      enddo
      close(334)
      
      inquire(unit = 335, opened = taAberta)
      if(taAberta)close(335)
      open(335, file=trim(adjustl(local))//trim(adjustl('lambC_vs_Y4Est.dat')), access = 'append', status='unknown')
            
      lamb = lambC
      
      call condicao_inicial(rede, alp, lamb, mu, 0.0_dp, tTotal, tRex)
      
      call sirs_estoc(rede, alp, lamb, mu, T_vs)
      write(335,*) lamb, Y4
      
      stop "Concluiu o calculo do IPR no limiar"
         
   else
      open(333, file=trim(adjustl(local))//trim(adjustl('lbd_vs_rhoEst.dat')), access='append', status='unknown')
      open(334, file=trim(adjustl(local))//trim(adjustl('lbd_vs_XiEst.dat')), access='append', status='unknown')
      open(335, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Y4Est.dat')), access='append', status='unknown')
      open(336, file=trim(adjustl(local))//trim(adjustl('lbd_vs_t_LS.dat')), access = 'append', status='unknown')
   endif
!#######################################################################
   
   do i2 = 1, nlamb
!#######################################################################   
      if(lamb >= lambdaf) exit
      call condicao_inicial(rede, alp, lamb, mu, 0.0_dp, tTotal, tRex)
      
      call sirs_estoc(rede, alp, lamb, mu, T_vs)
      write(333,*) lamb, rho_medioQS_global
      write(334,*) lamb, Xi_global
      write(335,*) lamb, Y4
      write(336,*) lamb, t_LS

!#######################################################################
      lamb = lamb + dlamb
!#######################################################################      
   enddo
   close(333); close(334); close(335); close(336)
!#######################################################################   

  contains

      subroutine criaRedeEClassificaClusters(this, tam_rede1, grau_min1, grau_max1, gama_exp1, seme1)
        type(grafo_PL_UCM) :: this
        integer :: tam_rede1
        integer :: grau_min1
        real(dp) :: grau_max1
        real(dp) :: gama_exp1
        integer :: seme1
!#######################################################################   
        call this%iniciaGrafo(tam_rede1)
!#######################################################################
        call this%inicia(grau_min1, grau_max1, gama_exp1, seme1)
!#######################################################################
        call this%liga(seme1, .False.) 
!#######################################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
!#######################################################################
      end subroutine
!#######################################################################

   subroutine calcula_Pk_E_Lambda0(this)
      class(grafo) :: this
!#######################################################################
      if(allocated(P_grau)) deallocate(P_grau)
      allocate(P_grau(this%degMin:this%degMax))
      P_grau = 0d0
      do i1 = 1, this%nodes
         P_grau(this%deg(i1)) = P_grau(this%deg(i1)) + 1.0d0
      enddo
      P_grau = P_grau/(1d0 * this%nodes)
      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//trim(adjustl(indice))//'.dat'
      open(800, file=trim(adjustl(arq_1)), status='unknown')
      qm = 0.0_dp; q2m = 0.0_dp
      do i1 = this%degMin, this%degMax
         if(P_grau(i1) == 0.0_dp) cycle
         write(800,*) i1, P_grau(i1)
         qm = qm + (1.0_dp * i1) * P_grau(i1)
         q2m = q2m + P_grau(i1) * ((1.0_dp * i1)**2.0_dp)
      enddo
      if ( lamb0 == 0.0_dp )then
         !call limiarQMF(rede, lamb0)
         !write(*,*) "Limiar QMF via subroutine limiarQMF: lambC = ", lamb0
         lamb0 = min((qm/q2m), ((1.0_dp * this%degMax)**(-0.5_dp)))
         write(*,*) "Limiar QMF a ser utilizado eh: ", lamb0
         lamb0 = lamb0 - 10.0d0 * dlamb
         write(*,*) "Lamb0 estimado via limiar QMF eh = ", lamb0
      endif
      nlamb = int((lambdaf-lamb0)/dlamb)
      write(*,*) "nlamb = ", nlamb      
      
      close(800)
      deallocate(P_grau)      
   end subroutine   
  
!=======================================================================
   subroutine kNN_e_clustering(this)
      class(grafo) :: this

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call calcula_k_nn(this, .True., 800, trim(adjustl(arq_1)))
      close(800)

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call clustering(this,.True., 800, trim(adjustl(arq_1)))
      close(800)
   end subroutine
!======================================================================= 
     
   subroutine acha_cutoff_rigido(kming, gam_p, N_p)
      integer :: kming
      real(dp) :: gam_p
      integer :: N_p
      real(dp), parameter :: tol = 5d-5
      real(dp) :: gminus, gmais
      real(dp) :: kminus, kmais
      real(dp) :: k_p
      real(dp) :: gp
      real(dp) :: Ap1
      integer :: kl1, kl2
      integer, parameter :: N_it = 10**3
 
      !#################################################################
      !   Inicio
      !#################################################################           
      kminus = 1d0 * kming
      kmais = 1.5d0 * kming * (1d0 * N_p)**(1d0/gam_p)
      
      if(allocated(Ap_list)) deallocate(Ap_list)
      allocate(Ap_list(int(kming):(int(kmais))))
      
      gminus = g_func(kminus, gam_p, N_p)
      
      if(gminus >= 0d0) stop "Precisamos de um valor de kminus para que gminus < 0"
       
      Ap1 = 0d0
      do kl1 = kming, int(kmais)
         Ap1 = Ap1 + (1d0 * kl1)**(-gam_p)
         Ap_list(kl1) = Ap1
      enddo

      gmais = g_func(kmais, gam_p, N_p)
   
      if(gmais <= 0d0) stop "Precisamos de um valor para kmais, tal que gmais > 0"
      !#################################################################
      !   Execucao
      !#################################################################
      kl1 = 1
      do while(kl1 <= N_it)
         k_p = kminus + (kmais - kminus)/2d0
         gp = g_func(k_p, gam_p, N_p)
         if((gp == 0d0).or.((kmais-kminus)/2d0 <= tol))then
            grau_max = k_p
            write(*,*) "Achou o grau max apos N = ", kl1, " iteracoes."
            exit
         endif
         kl1 = kl1 + 1
         if(gminus * gp > 0d0)then
            kminus = k_p
            gminus = gp
         else
            kmais = k_p
         endif
      enddo      
      !#################################################################
      !   Final
      !#################################################################
      deallocate(Ap_list)                   
   end subroutine
   
   function g_func(k_s, gam1, Nstr)
      real(dp) :: k_s
      real(dp) :: gam1
      integer :: Nstr
      real(dp) :: g_func
      real(dp) :: Ap
      !#################################################################
      !   Quando g_func for usada pela primeira vez,
      !   gbuffer = 0.0d0 e kmin = this%degMin,
      !   k_s recebe o valor que se quer aproximar de kc
      !   
      !#################################################################      
      !g_func = gbuffer
                  
      Ap = Ap_list(int(k_s))
      
      Ap = 1d0/Ap
      
      g_func = k_s - Ap**(1d0/gam1) * (1d0 * Nstr)**(1d0/gam1)
      
   end function 
   
end program
