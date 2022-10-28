!include ''
module sirs_estocastico
   use geraRede
   use mod_rndgen
   use mod_tools
   !use mod_tictoc
   implicit none
   
   integer, parameter :: dp = kind(0.0d0)
   !====================================================================    
        real(dp), allocatable :: Pn_QS_global(:), Pn_Rec_QS_global(:)
	real(dp):: rho_medioQS_global
        real(dp) :: rec_medioQS_global
	real(dp) :: rho2_medioQS_global
        real(dp) :: rec2_medioQS_global
	real(dp) :: dev_rho_medioQS_global 
	real(dp) :: Xi_global
        real(dp) :: Xi_rec_global
	real(dp) :: S_Schanon_global
        real(dp) :: S_rec_Shanon_global
	real(dp) :: t_LS
	real(dp) :: Y4
	real(dp), allocatable :: v1(:)
	integer :: ind_v1
	real(dp), parameter :: qtdade = 1.0d0
	integer :: stubs_infTotal
   integer(kind = 8) :: int_8
   integer :: ind_hub
   !====================================================================      
!   real(dp), allocatable :: I_i(:)
!   real(dp), allocatable :: R_i(:)
   !====================================================================
!  Tem tamanho this%list e guarda os estados 0, 1, 2 para cada sitio         
!   integer(kind=1), allocatable :: sigma(:)
   !====================================================================
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
!#######################################################################
   integer :: nInf
   integer :: nRec
   integer :: nSus   
   integer, allocatable :: inf_List(:)
   integer, allocatable :: rec_List(:)
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
!  ens_Inf_List(:,:) eh uma matriz na qual cada coluna representa
!  uma configuracao e cada linha representa sitios no estado infectado
!  naquela configuração.
!  Analogamente para a matriz ens_Rec_List(:,:).
!  Note, n_ens eh o numero maximo!
!  Em ens_nInf(:), um vetor, guardamos o numero de sitios infectados
!  naquela configuracao. A ens_nInf(i1), membro i1 da colecao,
!  corresponde uma confugracao ens_Inf_List(1:ens_nInf(i1), i1)
!  com ens_nInf(i1) sitios infectados.
!  O analogo ocorre com ens_Rec_List(:,:).
!####################################################################### 
   !==================================================================== 
   integer :: n_ens
   integer :: residuo
   integer, allocatable :: ens_Inf_List(:,:)
   integer, allocatable :: ens_Rec_List(:,:)
   integer, allocatable :: ens_nInf(:)
   integer, allocatable :: ens_nRec(:)
   !====================================================================
   real(dp) :: t, dt, t_relax, t_max
   !==================================================================== 
   integer :: nInfMax
   !character(len=500), parameter :: local = trim(adjustl("/home/jota/SIRS_Estocastico/Rst_Rede_Real/"))
   integer(kind=8) :: n_ite
   real(dp) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
   real(dp) :: t1_T, tf_T
   real(dp) :: del_t1, del_t2, del_t3, del_t4, del_t5   
   logical :: contaTempo   
   !==================================================================== 
   !type(tictoc) :: rel1, rel2, rel3, rel4, rel5, rel6, rel7, rel8, rel9, rel10, rel11
   !==================================================================== 
   contains

   !==================================================================== 
   subroutine aloca_listas_dinamicas(this)
      class(grafo), intent(in) :: this
      integer :: j1, j2
   !==================================================================== 
      if(allocated(Pn_QS_global)) deallocate(Pn_QS_global)
         allocate(Pn_QS_global(this%nodes))
   !====================================================================
      if(allocated(Pn_Rec_QS_global)) deallocate(Pn_Rec_QS_global)
         allocate(Pn_Rec_QS_global(0:this%nodes))    
   !==================================================================== 
      if(allocated(sigma)) deallocate(sigma)
         allocate(sigma(this%nodes))
   !==================================================================== 
      if(allocated(inf_List)) deallocate(inf_List)
         allocate(inf_List(this%nodes))
   !==================================================================== 
      if(allocated(rec_List)) deallocate(rec_List)
         allocate(rec_List(this%nodes))
   !==================================================================== 
         n_ens = 100
         residuo = int( qtdade * this%nodes)
   !==================================================================== 
   if(allocated(ens_Inf_List)) deallocate(ens_Inf_List)
      allocate(ens_Inf_List(residuo, n_ens))
   !==================================================================== 
   if(allocated(ens_nInf)) deallocate(ens_nInf)
      allocate(ens_nInf(n_ens))
   !==================================================================== 
   if(allocated(ens_Rec_List)) deallocate(ens_Rec_List)
      allocate(ens_Rec_List(residuo, n_ens))
   !==================================================================== 
   if(allocated(ens_nRec)) deallocate(ens_nRec)
      allocate(ens_nRec(n_ens))
   !==================================================================== 
      stubs_infTotal = 0
      do j1 = 1, this%nodes
         stubs_infTotal = stubs_infTotal + this%deg(j1)
      enddo
   !==================================================================== 
      if(allocated(v1)) deallocate(v1)
         allocate(v1(this%nodes))
      l_search_hub:do j1 = 1, this%nodes
         if( this%deg(j1) == this%degMax)then
            ind_hub = j1
            exit l_search_hub
         endif
      enddo l_search_hub
   !#######################################################################         
   end subroutine
   
   !====================================================================
   subroutine condicao_inicial(this, alp, lamb, mu, t_0, tMax, tRelax)
      !integer, parameter :: dp = kind(0.0d0)
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp, lamb, mu
      real(dp), intent(in) :: t_0, tMax, tRelax
      integer :: i1, l10
      integer :: j1
      integer(kind = 8) :: j12
      !=================================================================
      nInfMax = this%nodes
      nInf = nInfMax
      nRec = 0
      nSus = 0
      !=================================================================
      !  Vamos infectar todo mundo
      !=================================================================
      do i1 = 1, this%nodes
         if( lista_de_clusters(i1) /= i_comp_gigante )then
            sigma(i1) = 3
            cycle
         endif
         
         sigma(i1) = 1
         inf_List(i1) = i1
      enddo
      !=================================================================
      v1 = 0.0d0       
      !=================================================================
      !rec_List = 0      
      Pn_QS_global = 0.0_dp
      Pn_Rec_QS_global = 0.0_dp
      !=================================================================
      t = t_0
      t_max = tMax
      t_relax = tRelax 
      !=================================================================            
   end subroutine
   !====================================================================   
   subroutine sirs_estoc(this, alp, lamb, mu, T_vs)
      class(grafo), intent(in) :: this
      real(dp), intent(in) :: alp, lamb, mu
      logical, intent(in) :: T_vs
      real(dp) :: Tax_inf, Tax_rec, Tax_wan
      real(dp) :: Tax_total
      real(dp) :: prob_inf, prob_rec, prob_wan
      logical :: recup, recai, infec
      integer :: stubs_inf
      real(dp) :: prob
      type(rndgen) :: gen
      integer :: seed
      integer :: ultima_foto

      integer, parameter :: nfotos = 100
      integer :: i1, i2, i3, j1, j2, j3, k1, k2
      integer :: coluna
      integer(kind=8) :: j12
      real(dp) :: probInf, prob_dist_t
      real(dp), parameter :: p_deco = 0.02_dp
      real(dp) :: prob1, p_act, p_fill
      integer :: n_vezes
      integer :: sumPQS      
      character(len=100) :: dados
      character(len=10) :: alp_char, lamb_char
      character(len=50) :: rede_char
      integer :: fls_inf
      real(dp) :: tx_fls_inf
      integer :: cont_loop
      !real(dp), allocatable :: cont_loop_vec(:)
      character(len=9) :: tam_str
      !-----------------------------------------------------------------
      logical :: escr_dez_p

      real(dp) :: t_registrado

      integer :: um_p

      !-----------------------------------------------------------------
      !if(allocated(cont_loop_vec) ) deallocate(cont_loop_vec)
      !=================================================================
      !allocate(cont_loop_vec(0:this%nodes))
      !cont_loop_vec = 0.0d0      
      !=================================================================
      !call rel10%start()
      !call rel10%tic()
      write(alp_char,'(f8.2)') alp
      alp_char = trim(adjustl(alp_char))
      !=================================================================
      stubs_inf = stubs_infTotal
      ultima_foto = 0
      !=================================================================
      !contaTempo = .False.
      !=================================================================
      !if(contaTempo) call cpu_time(t1_T)
      !=================================================================

      um_p = int(t_max/100.0d0)

      t_registrado = 0.0d0

      write(*,'(A, 1X, I0, 1X, A)') 'Prog. = ', 0, '%'

      escr_dez_p = .True.

      !=================================================================
ld: do while(t <= t_max)

         !--------------------------------------------------------------
         ! Cast t to an 8 bytes variable.
         !--------------------------------------------------------------
         
         int_8 = t
         
         !--------------------------------------------------------------
         ! Check if progress was already registered.
         !--------------------------------------------------------------

         if( escr_dez_p .and. ( t > (1.0d0 * t_registrado + 5d0 * um_p) ) )then

            escr_dez_p = .False.

         elseif((.not. escr_dez_p) .and. ( (t > t_registrado + 10.0d0 * um_p) .and. (t < t_registrado + 15.0d0 * um_p) ))then

            t_registrado = t

            write(*,'(A, 1X, I0, 1X, A)') "Prog. = ", int(100.0d0 * (t/t_max)), "%"

            escr_dez_p = .True.

         endif
         !--------------------------------------------------------------
         !#######################################################################
         ! Calcula as taxas de eventos
         !#######################################################################        
         Tax_inf = 1.0_dp * stubs_inf * lamb
         Tax_rec = 1.0_dp * nInf * mu
         Tax_wan = 1.0_dp * nRec * alp     
         !-----------------------------------------------------------------------             
         
         Tax_total = Tax_inf + Tax_rec + Tax_wan
!#######################################################################
! Calcula as probabilidades de eventos acontecerem.
!#######################################################################
       prob_rec = 1.0_dp * Tax_rec/Tax_total
       prob_inf = 1.0_dp * Tax_inf/Tax_total
       prob_wan = 1.0_dp * Tax_wan/Tax_total
!#######################################################################
!   Calcula os passos de tempo, segundo o Processo de Poisson
!#######################################################################        
       prob = gen%rnd()
       dt = -1.0_dp * log(max(1e-12,prob))/Tax_total
       t = t + dt
!#######################################################################
! Pouco a pouco calcula a probabilidade quase-estacionária.
!#######################################################################
       if(t >= t_relax) Pn_QS_global(nInf) = Pn_QS_global(nInf) + dt
       if(t >= t_relax) Pn_Rec_QS_global(nRec) = Pn_Rec_QS_global(nRec) + dt
!#######################################################################
! Para t > 5.0, começamos a salvar configuracoes segundo
! pede o modelo QS
!#######################################################################
         !##############################################################
         !   Inspecionado 02/09   14:43
         !##############################################################
         if(t > 5_dp)then
            !===========================================================
            if(ultima_foto < nfotos)then
               prob1 = gen%rnd()
               p_fill = 1_dp * dt
               !========================================================
               if(prob1 <= p_fill)then                  
                     ultima_foto = ultima_foto + 1
                     ens_nInf(ultima_foto) = nInf
                     ens_nRec(ultima_foto) = nRec
                     !==================================================
                     do j1 = 1, nInf
                        ens_Inf_List(j1, ultima_foto) = inf_List(j1)
                     enddo
                     !==================================================
                     do j1 = 1, nRec                   
                        ens_Rec_List(j1, ultima_foto) = rec_List(j1)                                                                   
                     enddo                
                     !==================================================
!#######################################################################
               endif
!#######################################################################
            else
               prob1 = gen%rnd()
               p_act = 1.0_dp * p_deco * dt
!#######################################################################               					
               if(prob1 <= p_act)then
                     !==================================================               
                     k1 = gen%int(1, ultima_foto)
                     !==================================================
                     ens_nInf(k1) = nInf
                     !==================================================                     
                     do j1 = 1, nInf
                        ens_Inf_List(j1, k1) = inf_List(j1)
                     enddo
                     !==================================================                     
                     ens_nRec(k1) = nRec
                     !==================================================                     
                     do j1 = 1, nRec
                        ens_Rec_List(j1, k1) = rec_List(j1)                        
                     enddo
                     !==================================================                     
               endif
!#######################################################################
            endif
!#######################################################################
         endif
!#######################################################################
 ! Aqui acontecem as transicoes
!#######################################################################
       prob = gen%rnd()
!#######################################################################      
       if(prob <= prob_rec)then
            !###########################################################
            ! Aqui tah inspecionado jah.    02/09/21 14:07
            !###########################################################            
            if(t < 20.0_dp)then
               !########################################################
               ! Se t < 20.0, usa-se a condicao de contorno refletora
               !########################################################              
               if(nInf == 1)then
                  !#####################################################
                  ! O ultimo sitio infectado, que iria se recuperar.
                  !#####################################################
                  !=====================================================
                  j1 = inf_List(nInf) 
                  !=====================================================
                  stubs_inf = stubs_inf - this%deg(j1)
                  !=====================================================
                  sigma(j1) = 2
                  nRec = nRec + 1
                  rec_List(nRec) = j1                       
                  !=====================================================
                  !#####################################################
                  ! Sorteia numero aleatorio para escolher sitios
                  ! ou recuperados ou suscetiveis.
                  ! Fiz isso pq estava muito complicado de escolher
                  ! um sitio aleatoriamente e saber onde, na lista
                  ! de recuperados ele estaria.
                  !#####################################################
                  prob = gen%rnd()
                  !#####################################################
                  ! if(True): Se um sitio recuperado for selecionado,
                  ! escolhemos um recuperado aleatorio da lista.
                  ! else: um sitio da rede eh selecionado,
                  ! caso ele seja suscetivel, o aceitamos e prosseguimos
                  !#####################################################
                  if(prob <= (1.0_dp * nRec)/(1.0d0 * this%nodes))then                     
                     j3 = gen%int(1, nRec)
                     j1 = rec_List(j3)                     
                     rec_List(j3) = rec_List(nRec)                   
                     nRec = nRec - 1
                  else
                     lant:do
                       j1 = gen%int(1, this%nodes)
                       if(sigma(j1) == 0) exit lant
                     enddo lant
                  endif
                  !=====================================================
                  inf_List(nInf) = j1
                  sigma(j1) = 1
                  stubs_inf = this%deg(j1)
                  !=====================================================
                  cycle ld
               endif
            endif
            !###########################################################   
            ! * Um sitio infectado eh selecionado ao acaso.
            !###########################################################
            ! * O ultimo sitio da lista de infectados (inf_List)\
            !   eh colocado na posicao do sitio recem recuperado;
            ! * O numero de infectados eh reduzido em uma unidade;
            ! * O numero de stubs infectantes eh reduzido\
            !   em uma quantidade igual ao grau do sitio recem\
            !   recuperado.
            !###########################################################
            !###########################################################
            ! Aqui tah inspecionado jah.   02/09/21 14:07
            !###########################################################                        
            !===========================================================
            j3 = gen%int(1, nInf)
            j1 = inf_List(j3)
            !===========================================================
            sigma(j1) = 2
            nRec = nRec + 1
            rec_List(nRec) = j1            
            !===========================================================
            inf_List(j3) = inf_List(nInf)
            nInf = nInf - 1
            !===========================================================                        
            stubs_inf = stubs_inf - this%deg(j1)                    
            !==============================================================
            ! * Seu intervalo de atividade na dinamica eh incrementado.
            !==============================================================
            if(t >= t_relax) v1(j1) = v1(j1) + 1.0_dp/mu
            !==============================================================
            ! * Seu status eh mudado para 'recuperado' (#2);
            ! * O numero de recuperados eh incrementado;
            ! * O novo sitio recuperado eh adicionado aa lista re_List
            !==============================================================
            ! * Se o sitio recem recuperado eh uma folha,
            !   reduzimos em uma unidade o numero de folhas infectadas;
            ! * Testamos se a contagem estah sendo feita corretamente.
            !==============================================================

            !==============================================================
            ! # Se o estado absorvente eh alcancado, escolhemos\
            !  uma configuracao ja visitada;
            !==============================================================
            !   Aqui parece bem inspecionado ja tambem    02/09/21 14:10
            !==============================================================
                        
			   if(nInf == 0)then
               !========================================================
               if (stubs_inf /= 0) stop "Wrong counting in line 466."
               !========================================================
			      ! 1. Colocamos como *suscetiveis* (#0) todos os sitios\
               !    ainda recuperados;
               !========================================================
               do j3 = 1, nRec
                  sigma(rec_List(j3)) = 0
               enddo
               !########################################################
               ! 2. Escolhemos uma configuracao passada salva e\
               !    atribuimos a ela o indice *k1*;
               !========================================================
				   k1 = gen%int(1,ultima_foto)
               !========================================================
               ! 3. A esta configuracao corresponde o numero\
               !    de infectados *nInf* = *ens_nInf(k1)*;
               !========================================================
				   nInf = ens_nInf(k1)
               !========================================================
               ! 4. Atribuimos aa lista de infectados atual *inf_List*\
               !    cada um dos *nInf* sitios infectados salvos\
               !    na lista *ens_Inf_List(1:nInf, k1)*;
               ! 5. Mudamos o status de cada sitio atribuido para\
               !    *sigma(j1)* = 1;
               ! 6. Incrementamos o numero de stabs infectantes\
               !    em uma quantidade igual ao grau *this%deg(j1)*;
               ! 7. Se o sitio for uma folha, incrementamos o numero\
               !    de folhas infectantes em uma unidade\
               !    e testamos se a contagem estah correta.
               !========================================================         
				   do j3 = 1, nInf                                     
				      
                  j1 = ens_Inf_List(j3, k1)

                  inf_List(j3) = j1
				      
                  sigma(j1) = 1
				      
                  stubs_inf = stubs_inf + this%deg(j1) 
               
               enddo
               !========================================================
               ! 27/08/21
               !
               ! Aqui costumava aparecer indices negativos para\
               ! os sitios
               !
               !========================================================
				   ! 8. Atualizamos o numero de recuperados *nRec* atual\
				   !    para o numero de recuperados *ens_nRec(k1)*\
				   !    associado aa configuracao *k1*
               !========================================================
				    nRec = ens_nRec(k1)
				    !#######################################################
				    ! 9. Atribuimos aa lista de recuperados *rec_List*\
				    !    cada um dos *nRec* sitios recuperados contidos\
				    !    na lista *ens_Rec_List(1:nRec,k1)*;
				    ! 10. O seu status eh atualizado para *recuperado* (#2)\
				    !     na lista *sigma*
				    !=======================================================
				    ! *Por enquanto, eh aqui que estah aparecendo indice\
				    ! nulo na lista sigma. Onde sera a raiz disso?*
				    !=======================================================				
				   
				    !#######################################################
				    do j3 = 1, nRec  
				       j1 = ens_Rec_List(j3, k1)
				       rec_List(j3) = j1
				       sigma(j1) = 2
				    enddo
			    endif
         elseif(prob <= (prob_rec+prob_inf))then
               !########################################################
               ! Aqui realizamos o evento de infeccao.
               !########################################################
               ! Aqui tah acontecia muita rejeicao de folhas, por isso\
               ! o codigo costumava ser lento.
               !########################################################
               
               !########################################################
               ! Implementacao do Wesley para achar o lifespan
               !########################################################
               !call cpu_time(t1)
               !write(*,*) "t1 = ", t1
               !--------------------------------------------------------
               ! If the hub is infected, we consider the heterogeneity
               ! present in the infective degree of the nodes by means
               ! of the probability
               ! \Pi = \frac{this%deg(j1)}{\sum_i this%deg(i)}.
               ! But, if the hub is not infected, we know that all the
               ! infective nodes have the same degree of infectivity,
               ! so we choose the node homogeneously.
               !--------------------------------------------------------
               if(sigma(ind_hub) == 1)then
                  !-----------------------------------------------------
                  prob = gen%rnd()
                  if( prob <= 1.0d0 * (nInf - 1) * this%degMin/(1.0d0 * (this%degMax + (nInf - 1) * this%degMin ) ) )then
                     sel1: do
                        j3 = gen%int(1,nInf)
                        j1 = inf_List(j3)
                        if( j1 /= ind_hub ) exit sel1
                     enddo sel1
                  else
                     j1 = ind_hub
                  endif
                  !-----------------------------------------------------
               else
                  j3 = gen%int(1,nInf)
                  j1 = inf_List(j3)
               endif
               !========================================================              
               j12 = gen%int(this%aux(j1), this%aux(j1) + this%deg(j1) - 1)
               j2 = this%listAdj(j12)
               !========================================================               
               if(sigma(j2) == 0)then
                  !=====================================================
                  sigma(j2) = 1
                  stubs_inf = stubs_inf + this%deg(j2)                  
                  !=====================================================
                  nInf = nInf + 1
                  inf_List(nInf) = j2
                  !=====================================================
               endif
       else
          !=============================================================
          j3 = gen%int(1,nRec)
          j1 = rec_List(j3)
          !=============================================================
          sigma(j1) = 0
          !=============================================================
          rec_List(j3) = rec_List(nRec)
          nRec = nRec - 1
          !=============================================================
       endif
    enddo ld

   !====================================================================
   !if(contaTempo) call cpu_time(tf_T)
   !if(contaTempo) write(*,*) "Tempo total requerido foi ", tf_T - t1_T
   !====================================================================
   if(T_vs)then
      close(888)
   endif
   !====================================================================
	sumPQS = sum(Pn_QS_global)
   if( sumPQS <= 0.0d0) stop "Something scrued the quasi-stationary probability." 
   !====================================================================
	do i1 = 1, size(Pn_QS_global)			
		Pn_QS_global(i1) = 1.0_dp * Pn_QS_global(i1) / sumPQS
        Pn_Rec_QS_global(i1) = 1.0_dp * Pn_Rec_QS_global(i1) / sumPQS
	enddo
   !====================================================================
	rho_medioQS_global = 0.0_dp
	rho2_medioQS_global = 0.0_dp
	S_Schanon_global = 0.0_dp
   !====================================================================
    rec_medioQS_global = 0.0_dp
    rec2_medioQS_global = 0.0_dp
    S_rec_Shanon_global = 0.0_dp
   !====================================================================
	if(Pn_QS_global(1) > 0.0_dp)then
		t_LS = 1.0_dp /Pn_QS_global(1)
	else
		t_LS = 10.0d0**9.0d0
	endif
   !====================================================================
	do i1 = 1, nInfMax
           if(Pn_QS_global(i1) > 0.0_dp)then
              rho_medioQS_global = 1.0_dp * i1 * Pn_QS_global(i1) + rho_medioQS_global						! n medio
              rho2_medioQS_global =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_global(i1) + rho2_medioQS_global
              S_Schanon_global = S_Schanon_global - Pn_QS_global(i1) * log(Pn_QS_global(i1))
           endif

           if(Pn_Rec_QS_global(i1) > 0.0_dp)then
              rec_medioQS_global = 1.0_dp * i1 * Pn_Rec_QS_global(i1) + rec_medioQS_global
              rec2_medioQS_global =  ((1.0_dp *i1) ** 2.0_dp) * Pn_Rec_QS_global(i1) + rec2_medioQS_global
              S_Rec_Shanon_global = S_Rec_Shanon_global - Pn_Rec_QS_global(i1) * log(Pn_Rec_QS_global(i1))
           endif
	enddo

    rho_medioQS_global = 1.0_dp * rho_medioQS_global /(1.0d0 * this%nodes)
    
    rho2_medioQS_global = 1.0_dp * rho2_medioQS_global/ ((1.0_dp * this%nodes)**2.0_dp)

    dev_rho_medioQS_global = (rho2_medioQS_global - rho_medioQS_global**2.0_dp)**0.5_dp

    Xi_global = 1.0_dp * this%nodes * (rho2_medioQS_global - (rho_medioQS_global**2.0_dp))/rho_medioQS_global


    rec_medioQS_global = 1.0_dp * rec_medioQS_global / (1.0d0*this%nodes)

    rec2_medioQS_global = 1.0_dp * rec2_medioQS_global/ ((1.0_dp * this%nodes)**2.0_dp) 

    Xi_Rec_global = 1.0_dp * this%nodes * (rec2_medioQS_global - (rec_medioQS_global**2.0_dp))/rec_medioQS_global


	v1 = v1/(1.0_dp * (t - t_relax))
	
	v1 = v1/(sum(v1**2.0_dp))**0.5_dp

	Y4 = sum(v1**4.0_dp)
    
        !###############
        !call rel9%toc()
        !###############
    
    write(lamb_char, '(F8.5)') lamb
    
    write(tam_str, '(I0)')this%nodes
    
    !open(998, file = 'N_rejeicoes_sel_infeccioso_lamb_'//trim(adjustl(lamb_char))//'_tam_'//trim(adjustl(tam_str))//'.dat', status = 'unknown')
    
    !cont_loop_vec = cont_loop_vec/(1.0d0 * sum(cont_loop_vec) )
    
    !do j1 = 0, size(cont_loop_vec)
       !if(cont_loop_vec(j1) > 0.0d0) write(998, *) j1, cont_loop_vec(j1)
    !enddo
    !deallocate(cont_loop_vec)
    !close(998)
   end subroutine
!#######################################################################
   
end module


program main
   use geraRede
   use sirs_estocastico
   use mod_rndgen
   use mod_tools
   implicit none
   !====================================================================  
   type(grafoRRN_Plus_Star) :: rede
   !type(grafo_Erdos_Renyi) :: rede
   !====================================================================   
   real(dp) :: dlamb
   real(dp) :: lamb, lambdaf
   integer :: ind_lamb, n_lamb 
   !====================================================================      
   !Mudar aqui soh
   real(dp) :: lamb0 
   real(dp) :: alp
   !====================================================================
   real(dp), parameter :: mu = 1.0_dp
   !alp = 100 para testar sirs => sis
   !====================================================================
   type(rndgen) :: gen
   integer :: seed
   !====================================================================
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
   !====================================================================
   integer, allocatable :: tam_m1(:)
   !====================================================================   
   integer :: i1, i2, i3, j1, j2, j3, m10
   logical :: T_vs
   !====================================================================
   integer :: sumDeg2
   !====================================================================   
   real(dp) :: S_Schanon_C
   real(dp) :: Xi_max
   real(dp) :: Xi
   real(dp) :: lambC
   real(dp) :: t_LS_C
   !====================================================================
   character(len=300) :: cwd, resultados, tipoCorte
   character(len=1000) :: local
   character(len=1000) :: nomeArquivo
   character(len=20) :: buffer
   !==================================================================== 
   character(len=10) :: alp_char2  
   character(len=500) :: t_vs_Im
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo_rede_real
   character(len=7) :: tam_char
   character(len=5) :: grau_char, grauStar_char
   character(len=5) :: indice
   character(len=500) ::arq_1   
   character(len=3) :: char_ind_lamb
   character(len=10) :: lamb_char
   real(dp) :: Y4C
   !====================================================================
   integer :: nargus
   real(dp) :: divisor   
   real(dp) :: tMax1, tRelax1
   integer :: contador
   logical :: naoEscreveu
   integer :: interpola
   integer :: status_io
   integer :: grauRRN, grauStar
   character(len=1000) :: fi_char
   integer :: Eh_Limiar
!=======================================================================   
   resultados = trim(adjustl('Rst_Grafo_RRN_Mais_Estrela'))
   call system('mkdir -p '//resultados)
!=======================================================================   
   local = trim(adjustl(resultados))//"/"
!=======================================================================   
   call entradaArgumentos()
!=======================================================================
   local = trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//"/"

   call system('mkdir -p '//trim(adjustl(local)) )
   
   write(alp_char2, '(f5.1)') alp   
   alp_char2 = trim(adjustl(alp_char2))

   local = trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//"/"   

   call system('mkdir -p '//trim(adjustl(local)) )
   
!=======================================================================
   seed  = 967891968
!=======================================================================
   if( ind_lamb < -1)then
      stop "Valor de ind_lamb nao permite determinar lambda"
   elseif( ind_lamb == -1)then
      if(interpola /= -1) open(335, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', status='unknown')
   endif
!=======================================================================

!#######################################################################
! Eu setto Xi_max aqui e procuro nos arquivos. Se eu encontro um Xi_max,
! mesmo que ele nao seja o pico, ele serah util no trecho que calcula
! o Xi_max real.
!#######################################################################

!=======================================================================
   Xi_max = 0.0d0
   contador = 0
!=======================================================================
   if( interpola > 0)then      
      do
         read(335,*, iostat = status_io) lamb, Xi
         if(status_io /= 0)exit
         if( Xi > Xi_max)then
            lamb0 = lamb            
            Xi_max = Xi
            contador = 0
         else
            contador = contador + 1
         endif
      enddo
      close(335)  
   elseif( interpola == 0 )then
      !lamb0 = 0.0d0
      do
         read(335,*, iostat = status_io) lamb, Xi
         if(status_io /= 0) exit
         !lamb0 = 0.0d0
         !--------------------------------------------------------------
         if(Xi > Xi_max)then
            Xi_max = Xi
            contador = 0
         else   
            contador = contador + 1
         endif
         !--------------------------------------------------------------
         if( lamb > lamb0 ) lamb0 = lamb
         !--------------------------------------------------------------         
      enddo
      !stop "Leitura bem sucedida"
      lamb0 = lamb0 + dlamb
      close(335)
   endif
   !====================================================================
   if( interpola == 1)then
      !-----------------------------------------------------------------
      ! 3.5 passos para tras
      !-----------------------------------------------------------------
      lamb0 = lamb0 - 7.0d0 * dlamb/2.0d0 
   elseif( interpola == 2)then
      !-----------------------------------------------------------------
      ! 3.25 passos para tras
      !-----------------------------------------------------------------   
      lamb0 = lamb0 - (13.0d0/4.0d0) * dlamb
   elseif( interpola == 3)then
      !-----------------------------------------------------------------
      ! 3.75 passos para tras
      !-----------------------------------------------------------------   
      lamb0 = lamb0 - (15.0d0/4.0d0) * dlamb
   elseif( interpola == 4)then
      !-----------------------------------------------------------------
      ! 3.33 passos para tras
      !-----------------------------------------------------------------      
      lamb0 = lamb0 - (10.0d0/3.0d0) * dlamb
   elseif( interpola == 5)then
      !-----------------------------------------------------------------
      ! 3.67 passos para tras
      !-----------------------------------------------------------------         
      lamb0 = lamb0 - (11.0d0/3.0d0) * dlamb                    
   endif
   !====================================================================   
   !lamb = lamb0
   !====================================================================
   !open(335, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi'))//'_tam_'//trim(adjustl(tam_char))//'degStar_'//trim(adjustl(grauStar_char))//'degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', status='unknown')
   
   if(ind_lamb == (-1) )then
      write(*,*) "ind_lamb = ", ind_lamb
      open(333, file=trim(adjustl(local))//trim(adjustl('k_vs_t_LS'))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(334, file=trim(adjustl(local))//trim(adjustl('lbd_vs_rho'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(335, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(336, file=trim(adjustl(local))//trim(adjustl('lbd_vs_S_Shanon'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(337, file=trim(adjustl(local))//trim(adjustl('k_vs_S_ShanonC'))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(338, file=trim(adjustl(local))//trim(adjustl('k_vs_lbdC'))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(339, file=trim(adjustl(local))//trim(adjustl('k_vs_Y4C'))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(340, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Y4'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(341, file=trim(adjustl(local))//trim(adjustl('k_vs_Xi'))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')   
      !====================================================================
      open(342, file=trim(adjustl(local))//trim(adjustl('lbd_vs_rec'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(343, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi_rec'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      open(344, file=trim(adjustl(local))//trim(adjustl('lbd_vs_S_rec_Shanon'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'.dat', access = 'append', status='unknown')
      !-----------------------------------------------------------------
      lamb = lamb0      
      lambdaf = lamb0 + 500.0d0 * dlamb
      !-----------------------------------------------------------------
   elseif( ind_lamb >= 0)then
      write(*,*) "Ind_lamb = ", ind_lamb
      !-----------------------------------------------------------------
      lamb = lamb0 + 1.0d0 * ind_lamb * dlamb
      lambdaf = lamb
      !-----------------------------------------------------------------
      write(lamb_char, '(F8.5)') lamb
      !-----------------------------------------------------------------
      open(333, file=trim(adjustl(local))//trim(adjustl('k_vs_t_LS'))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(334, file=trim(adjustl(local))//trim(adjustl('lbd_vs_rho'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(335, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(336, file=trim(adjustl(local))//trim(adjustl('lbd_vs_S_Shanon'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(337, file=trim(adjustl(local))//trim(adjustl('k_vs_S_ShanonC'))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(338, file=trim(adjustl(local))//trim(adjustl('k_vs_lbdC'))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(339, file=trim(adjustl(local))//trim(adjustl('k_vs_Y4C'))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', access = 'append', status='unknown')
      open(340, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Y4'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', status='unknown')
      open(341, file=trim(adjustl(local))//trim(adjustl('k_vs_Xi'))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', status='unknown')
      !====================================================================
      open(342, file=trim(adjustl(local))//trim(adjustl('lbd_vs_rec'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', status='unknown')
      open(343, file=trim(adjustl(local))//trim(adjustl('lbd_vs_Xi_rec'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', status='unknown')
      open(344, file=trim(adjustl(local))//trim(adjustl('lbd_vs_S_rec_Shanon'))//'_tam_'//trim(adjustl(tam_char))//'_degStar_'//trim(adjustl(grauStar_char))//'_degRRN_'//trim(adjustl(grau_char))//'_alp_'//trim(adjustl(alp_char2))//'_lamb_'//trim(adjustl(lamb_char))//'.dat', status='unknown')   
   endif
   !====================================================================
   write(*,*) "O valor de alp eh ", alp
   write(*,*) "O valor de mu eh ", mu
   write(*,*) "O valor de lambda0 eh ", lamb
   !====================================================================
!#######################################################################
!				Inicia grafo
!#######################################################################
   !====================================================================
   !call rede%liga_estrela(tam_rede)
   call rede%ligaRRNStar(tam_rede, grauRRN, grauStar)
   !call rede%liga_ER(1000, 5.0d0)
   !====================================================================
   write(*,*) ""
   write(*,*) "O tamanho do grafo estrela eh", rede%nodes
   write(*,*)"O grau minimo da rede eh: ", rede%degMin
   write(*,*) ""
   write(*,*)"O grau maximo da rede eh: ", rede%degMax
   write(*,*) ""
   write(*,*)"O grau medio da rede eh: ", rede%degMean
   write(*,*) ""
   !====================================================================
   write(*,*) "Gerou a rede"
   write(*,*) ""
   !====================================================================
   call Dis_grau(rede)
   call criaRedeEClassificaClusters(rede)
   !==================================================================== 
   call aloca_listas_dinamicas(rede)
   !====================================================================
   !####################################################################
   ! Eu zero meu contador inicialmente na leitura de arquivos lah
   ! em cima, toda vez que um novo Xi_max aparece.
   ! Se nenhum Xi_max novo aparece, incremento o contador.
   ! De igual maneira eu zero meu contador na rotina abaixo quando
   ! acho um novo candidato a Xi_max
   ! se ele persistir por 5 vezes, ele leva o titulo.
   ! De inicio, coloco aqui um valor que nao permitira escrita caso
   ! jah haja um Xi_max.
   !####################################################################
   ! contador = 0
   !####################################################################   
   naoEscreveu = .True.
   
   ! Xi_max = 0.0d0. Mas ele eh settado lah em cima.
   
   do while(lamb <= lambdaf)
      call condicao_inicial(rede, alp, lamb, mu, 0.0_dp, tMax1, tRelax1)
      !=================================================================      
      call sirs_estoc(rede, alp, lamb, mu, T_vs)
      write(*,*) "-----------------------------------------------------"
      write(*,*) "Concluiu um ponto"
      if( (ind_lamb == 0) .and. (Eh_Limiar == 1) )then
         !--------------------------------------------------------------
         if( rede%nodes == 10**7 )then
            !-----------------------------------------------------------
            call calcula_distancias_ao_hub(rede)
            !-----------------------------------------------------------      
            fi_char = trim(adjustl(local))//trim(adjustl('fi_Sim_Stoch_N_10_a_7_RRN_deg_6_Plus_Star_998_alp_'//trim(adjustl(alp_char2))))
            !-----------------------------------------------------------
            open(999, file = trim(adjustl(fi_char))//'.dat', status = 'unknown')
            do i1 = 1, rede%nodes
               write(999,*) rede%deg(i1), lista_distancias(i1), (v1(i1))**2.0d0
               !write(*,*) "i1 = ", i1, "fi = ", v1(i1)
            enddo
            !-----------------------------------------------------------
            deallocate(lista_distancias)
            deallocate(v1)
            !-----------------------------------------------------------
            close(999)
            !-----------------------------------------------------------
            call system('tar -czf '//trim(adjustl(fi_char))//'.tar.gz '//trim(adjustl(fi_char))//'.dat' )
            call system('rm '//trim(adjustl(fi_char))//'.dat' )
            stop "Calculou apenas uma grandeza critica. Venha e me edite se quiser mais."
            !-----------------------------------------------------------      
         endif
      !-----------------------------------------------------------------      
      endif
      !=================================================================
      write(334,*) lamb, rho_medioQS_global
      !=================================================================
      write(335,*) lamb, Xi_global
      !=================================================================
      write(336,*) lamb, S_Schanon_global
      write(340,*) lamb, Y4
      !=================================================================  
      write(342,*) lamb, rec_medioQS_global
      write(343,*) lamb, Xi_rec_global
      write(344,*) lamb, S_rec_Shanon_global      
      !=================================================================            
      if( Xi_global > Xi_max )then
         Xi_max = Xi_global
         lambC = lamb
         S_Schanon_C = S_Schanon_global
         t_LS_C = t_LS
         Y4C = Y4
         contador = 0
      else
         contador = contador + 1
         if( (naoEscreveu) .and. (contador == 5) )then
            !===========================================================
            write(333,*) rede%degMax, t_LS_C
            write(337,*) rede%degMax, S_Schanon_C
            write(338,*) rede%degMax, lambC
            write(339,*) rede%degMax, Y4C
            write(341,*) rede%degMax, Xi_max
            !===========================================================
            naoEscreveu = .False.
            !===========================================================
            close(333)
            close(337)
            close(338)
            close(339)                         
            !===========================================================            
         endif
      endif
      !=================================================================      
      lamb = lamb + dlamb
      !=================================================================
   enddo   
   !====================================================================
   close(334)
   close(335)
   close(336)
   close(340)
   close(341)
   close(342)
   close(343)
   close(344)   
   !====================================================================  
!#######################################################################   
      !Processa os dados
      !Escreve os arquivos
!#######################################################################

contains

      subroutine entradaArgumentos()
         nargus = iargc()
         
         if(nargus == 12)then
            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede
            ! ${tam_rede}
            !#############################
            ! Lambda0
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauRRN
            write(grau_char, '(I0)') grauRRN
            ! ${grau_rrn}
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0 
            !  ${lamb0}           
            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
            ! ${divisor}
            !#############################
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor
            !#############################
            ! Alfa
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
            ! ${alp}
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) n_lamb
            ! ${n_lamb}
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tMax1
            ! ${tTotal}
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tRelax1
            !  ${tRelax}
            !#############################
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) interpola
            ! ${interpola}
            !#############################
            call getarg(10, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grauStar
            write(grauStar_char, '(I0)') grauStar
            ! ${grauStar}
            !#############################
            call getarg(11, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_lamb
            
            if(ind_lamb >= 0)then
               write(char_ind_lamb, '(I0)') ind_lamb
            else
               char_ind_lamb = ''
            endif

            call getarg(12, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) Eh_Limiar           
         else
            stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
         endif

         !==============================================================
         !  Sobrescrevo o lamb0 para 1/sqrt(tam_rede)
         !==============================================================
         if(interpola /= -1) lamb0 = 1.0d0/(1.0d0 * tam_rede)**0.5d0

         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         if(tam_rede == 10)then
            tam_char = '10'
         elseif(tam_rede == 100)then
            tam_char = '100'
         elseif(tam_rede == 1000)then
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
      end subroutine
!=======================================================================
      subroutine Dis_grau(this)
          !-------------------------------------------------------------
          class(grafo) :: this
          real(dp) :: P_k(0:this%nodes)
          integer :: i1 
          !-------------------------------------------------------------          
          P_k = 0.0d0
          
          do i1 = 1, this%nodes
             P_k(this%deg(i1)) = P_k(this%deg(i1)) + 1.0d0
          enddo
          
          P_k = P_k/(1.0d0 * this%nodes)
          
          do i1 = 1, this%nodes
             if(P_k(i1) > 0.0d0) write(*,*) i1, P_k(i1)
          enddo
                    
      end subroutine


      subroutine criaRedeEClassificaClusters(this)
        class(grafo) :: this
        !###############################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
        !###############################################################
      end subroutine
end program
