!include 'mod_rndgen_multiple.f90'
! include 'mod_tictoc.f90'
!#######################################################################
!	Modulo que gera a rede
!#######################################################################

module	geraRede
	use mod_rndgen
	implicit none

!#######################################################################
!	Parametros
!#######################################################################
	
	integer, private, parameter :: dp = kind(0.0d0)

!#######################################################################
!	Variaveis auxiliares
!#######################################################################	
	integer :: ki, k1, k2, N, kMin2
	real(dp) :: kMax2
	real(dp) :: p1, gama2, probTube, probFolhas
	integer :: id_amostra, seed_Ensemble
	integer, allocatable :: grau2(:)
	real(dp), allocatable :: pok(:)
	logical :: ressorteia
	integer :: sorteado
    character :: sorteado_char
!#######################################################################
!	Objeto tipo grafo
!#######################################################################

	type grafo
		integer :: nodes, degMin, degMax
		integer (kind=8) :: edge, sumDeg
		integer, allocatable:: deg(:), listAdj(:), aux(:), matriz(:,:)
		real(dp) :: degMean, degDev
		contains
			procedure :: iniciaGrafo
			procedure :: liga => ligaUCM
			procedure :: RedeReal => ligaRedeReal
	end type

    type, extends(grafo) :: grafo_Erdos_Renyi
       real(dp) :: prob_ER
       contains
          procedure :: liga_ER => liga_Erdos_Renyi
    end type
!#######################################################################
!	Objeto tipo grafo Lei de Potencia (Power Law) atraves do modelo UCM
!#######################################################################
	
	type, extends(grafo) :: grafo_estrela
	   contains
	      procedure :: liga_estrela => liga_grafo_estrela
	end type 

	type, extends(grafo) :: grafo_roda
	   contains
	   procedure :: liga_roda => liga_grafo_roda
	   procedure :: liga_roda_atacada => liga_grafo_roda_sob_ataque
	end type 

	
	type,extends(grafo) :: grafo_PL_UCM
		real(dp) :: gamm
		contains
			procedure ::	inicia => inicia_PL_UCM
	end type

!#######################################################################
!	Objeto tipo grafo Lei de Potencia (Power Law) atraves do modelo UCM
!	com controle de tubos
!#######################################################################
	
	type,extends(grafo_PL_UCM) :: grafo_PL_UCM_Tubes
		real(dp) :: prob_tube
		contains
			procedure ::	inicia_Tubes => inicia_PL_UCM_Tubes
	end type

!#######################################################################

	type,extends(grafo_PL_UCM) :: grafo_PL_UCM_Folhas
		real(dp) :: prob_folhas
		contains
			procedure ::	inicia_Folhas => inicia_PL_UCM_Folhas
	end type

!#######################################################################
!
!#######################################################################	
	
	type,extends(grafo) :: grafoRRN
	    integer :: degRRN
		contains
			procedure ::	iniciaRRN => iniciaGrafoRRN
	end type
	
	type,extends(grafo) :: grafoRRN_Plus_Star
	    integer :: degRRN, degStar
		contains
			procedure :: ligaRRNStar => ligaGrafoRRN_Plus_Star
	end type

!#######################################################################
!	Tive que por esta lista aqui porque este modulo eh comum
!	a outros dois modulos
!#######################################################################	
	
	integer, allocatable :: lista_de_clusters_din(:)
	integer(kind=1), allocatable :: sigma(:)
	
	real(dp) :: plus_xN, f_Tubos, alpha
	
	integer :: PLA_int
	logical :: PLA
	
	!###################################################################
	!		Subrotina que inicializa o grafo
	!###################################################################
	
	contains
	
	    subroutine liga_Erdos_Renyi(this, net_size, grau_medio)
	       !------------------------------------------------------------
	       class(grafo_Erdos_Renyi) :: this
	       integer :: net_size
	       real(dp) :: grau_medio
	       real(dp) :: prob_Erdos_Renyi
	       integer :: i1, i2
	       type(rndgen) :: ger
	       integer :: semente
	       real(dp) :: prob
	       !------------------------------------------------------------	       
	       prob_Erdos_Renyi = 1.0d0 * grau_medio/( 1.0d0 * net_size - 1.0d0 )
	       this%nodes = net_size
	       !------------------------------------------------------------
	       if(allocated(this%deg)) deallocate(this%deg)
	       allocate(this%deg(this%nodes))
	       !------------------------------------------------------------
	       if(allocated(this%aux) ) deallocate(this%aux)
	       allocate(this%aux(this%nodes))
	       !------------------------------------------------------------
	       semente = 98674890
	       call ger%init(semente)
	       this%deg = 0
	       !------------------------------------------------------------       
	       do i1 = 1, this%nodes
	          do i2 = i1 + 1, this%nodes
	             prob = ger%rnd()
	             if( prob <= prob_Erdos_Renyi )then
	                this%deg(i1) = this%deg(i1) + 1
	                this%deg(i2) = this%deg(i2) + 1
	             endif
	          enddo
	       enddo
	       !------------------------------------------------------------
	       this%aux(1) = 1
	       do i1 = 2, this%nodes
	          this%aux(i1) = this%aux(i1 - 1) + this%deg(i1 - 1)
	       enddo
	       !------------------------------------------------------------
	       this%sumDeg = sum(this%deg)
	       if(allocated(this%listAdj)) deallocate(this%listAdj)
	       allocate(this%listAdj(this%sumDeg))
	       !------------------------------------------------------------
	       semente = 98674890
	       call ger%init(semente)
	       this%deg = 0
	       !------------------------------------------------------------       
	       do i1 = 1, this%nodes
	          do i2 = i1 + 1, this%nodes
	             prob = ger%rnd()
	             if( prob <= prob_Erdos_Renyi )then
	                this%listAdj( this%aux(i1) + this%deg(i1) ) = i2
	                this%deg(i1) = this%deg(i1) + 1
	                this%listAdj( this%aux(i2) + this%deg(i2) ) = i1
	                this%deg(i2) = this%deg(i2) + 1
	             endif
	          enddo
	       enddo
	       !------------------------------------------------------------       
	    end subroutine
	
	
	    subroutine ligaGrafoRRN_Plus_Star(this, net_size, degRRN, degStar)
	       class(grafoRRN_Plus_Star) :: this
	       integer, intent(in) :: net_size
	       integer, intent(in) :: degRRN
	       integer :: degStar
	       integer :: i1, j1
	       integer :: stubs_disp
	       integer :: u_stub
	       integer(kind = 8) :: i2, i3, i4, stub1, stub2
	       integer :: site1, site2
	       integer :: semente
	       integer, allocatable :: degAux(:)
	       integer(kind = 8), allocatable :: lista_stubs(:)
	       type(rndgen) :: ger
	       
	       if( (net_size - 1) < degStar ) stop "Tamanho da rede nao permite conectar RRN ao hub. n_sitios_RRN < degStar."
	       !------------------------------------------------------------
	       ! Quando conectamos sitios rrn's ao hub, sobram 
	       ! degStar sitios com degRRN - 1 stubs disponiveis e !
	       ! (net_size - 1 - degStar) sitios com degRRN.
	       ! Entao, sobram 
	       ! degStar * (degRRN - 1) + (net_size - 1 - degStar) * degRRN
	       !------------------------------------------------------------
	       ! degStar * degRRB - degStar + net_size * degRRN - degRRN - degStar * degRRN
	       !------------------------------------------------------------
	       ! = net_size * degRRN - degStar - degRRN
	       ! = (net_size - 1) * degRRN - degStar
	       !------------------------------------------------------------
	       ! Define atributos iniciais da rede
	       !------------------------------------------------------------
	       this%nodes = net_size
	       !------------------------------------------------------------
	       ! Cria a lista de graus
	       if( allocated(this%deg) ) deallocate(this%deg)
	       allocate( this%deg(net_size) )
	       !------------------------------------------------------------
	       ! Se a soma de stubs dah impar, soma mais um stub ao hub.
	       !------------------------------------------------------------
	       stubs_disp = (net_size - 1) * degRRN - degStar
	       
	       this%degMax = degStar
	       if( mod(stubs_disp, 2) > 0 )then
	          this%degMax = this%degMax - 1
	       endif
	       !------------------------------------------------------------
	       ! Preenche a lista de graus
	       this%deg(1) = this%degMax
	       do i1 = 2, this%nodes
	          this%deg(i1) = degRRN
	       enddo
	       this%sumDeg = sum(this%deg)
	       !------------------------------------------------------------
	       if( allocated(degAux) ) deallocate( degAux )
	       allocate( degAux(this%nodes) )
           !------------------------------------------------------------
           degAux = 0
	       !------------------------------------------------------------
	       ! Cria a lista auxiliar
	       if( allocated(this%aux) ) deallocate(this%aux)
	       allocate( this%aux(net_size) )
	       !------------------------------------------------------------
	       ! Preenche a lista auxiliars
	       this%aux(1) = 1
	       do i1 = 2, this%nodes
	          this%aux(i1) = this%aux(i1 - 1) + this%deg(i1 - 1)
	       enddo
	       !------------------------------------------------------------
	       ! Cria a lista de adjacencias
	       if( allocated(this%listAdj) ) deallocate(this%listAdj)
	       allocate( this%listAdj(this%sumDeg) )
	       !------------------------------------------------------------
	       ! Conectamos o hub aos this%deg(1) primeiros sitios RRN.
	       !------------------------------------------------------------
	       do i1 = 2, this%deg(1) + 1
	          this%listAdj( this%aux(1) + degAux(1) ) = i1
	          degAux(1) = degAux(1) + 1
	          !---------------------------------------------------------
	          this%listAdj( this%aux(i1) + degAux(i1) ) = 1
	          degAux(i1) = degAux(i1) + 1
	          !---------------------------------------------------------
	       enddo
	       !------------------------------------------------------------
	       ! Cria lista de stubs
	       if( allocated(lista_stubs) ) deallocate( lista_stubs )
	       allocate( lista_stubs(this%sumDeg - 2 * this%deg(1) ) )
	       !------------------------------------------------------------
	       ! Preenche a lista de stubs com o restante dos stubs dos
	       ! sitios RRN.
	       !------------------------------------------------------------
	       u_stub = 0
	       do i1 = 2, this%nodes
	          j1 = degAux(i1)
	          do while( j1 < this%deg(i1) )
	             u_stub = u_stub + 1
	             lista_stubs(u_stub) = i1
	             j1 = j1 + 1
	          enddo
	       enddo
	       !------------------------------------------------------------
	       ! Aqui nos selecionamos os stubs ao acaso e os conectamos
	       !------------------------------------------------------------
	       semente = 98736476
	       call ger%init(semente)
	       !------------------------------------------------------------
lwhile:    do while( u_stub > 0)
	          stub1 = ger%int(1, u_stub)
	          site1 = lista_stubs(stub1)
              !---------------------------------------------------------
lstb2:        do
                 stub2 = ger%int(1, u_stub)
                 site2 = lista_stubs(stub2)
                 !if(site1 /= site2) exit lstb2
                 if( site1 /= site2 )then
                    do i1 = 1, degAux(site1)
                       if( this%listAdj(this%aux(site1) + i1 - 1) == site2 ) cycle lstb2
                    enddo                 
                    exit lstb2
                 endif
              enddo lstb2
              !---------------------------------------------------------
              !do i1 = 1, degAux(site1)
                 !if( this%listAdj(this%aux(site1) + i1 - 1) == site2 ) cycle lwhile
              !enddo
              !---------------------------------------------------------
              this%listAdj( this%aux(site1) + degAux(site1) ) = site2
              degAux(site1) = degAux(site1) + 1
              !---------------------------------------------------------
              this%listAdj(this%aux(site2) + degAux(site2)) = site1
              degAux(site2) = degAux(site2) + 1
              !---------------------------------------------------------
              if( stub1 > stub2)then
                 lista_stubs(stub1) = lista_stubs(u_stub)
                 u_stub = u_stub - 1
                 lista_stubs(stub2) = lista_stubs(u_stub)
                 u_stub = u_stub - 1
              else
                 lista_stubs(stub2) = lista_stubs(u_stub)
                 u_stub = u_stub - 1
                 lista_stubs(stub1) = lista_stubs(u_stub)
                 u_stub = u_stub - 1                 
              endif
	       !------------------------------------------------------------
	       enddo lwhile
	       
	       this%sumDeg = sum(this%deg)
	       this%degMin = minval(this%deg)
	       this%degMax = maxval(this%deg)
	       this%degMean = 1.0d0 * sum(this%deg)/(1.0d0 * this%nodes)
	       !------------------------------------------------------------
	       deallocate(lista_stubs, degAux)
	       !------------------------------------------------------------
	    end subroutine
!-----------------------------------------------------------------------
        subroutine liga_grafo_roda(this, grau_hub)
           class(grafo_roda) :: this
           integer, intent(in) :: grau_hub
           integer :: e1, e2
           
           this%nodes = grau_hub + 1
           
           if(allocated(this%deg)) deallocate(this%deg)
              allocate(this%deg(this%nodes))
           
           this%deg(1) = grau_hub
           
           do e1 = 2, this%nodes
              this%deg(e1) = 3
           enddo
           this%degMin = 3
           this%degMax = grau_hub
           this%sumDeg = 4 * grau_hub
           this%degMean = 1.0_dp * this%sumDeg/(1.0_dp * this%nodes)
           
           if(allocated(this%aux)) deallocate(this%aux)
              allocate(this%aux(this%nodes))           
           
           this%aux(1) = 1
           
           do e1 = 2, this%nodes
              this%aux(e1) = this%aux(e1-1) + this%deg(e1-1)
           enddo

           if(allocated(this%listAdj)) deallocate(this%listAdj)
              allocate(this%listAdj(this%sumDeg))           
           !------------------------------------------------------------
           ! Para o hub
           !------------------------------------------------------------
           do e1 = 1, this%degMax
              this%listAdj(e1) = e1 + 1
           enddo
           !------------------------------------------------------------
           ! /Para o hub
           !------------------------------------------------------------
           !------------------------------------------------------------
           ! Para os raios da roda
           !------------------------------------------------------------           
           this%listAdj(this%aux(2)) = 1
           this%listAdj(this%aux(2) + 1) = 3
           this%listAdj(this%aux(2)+2) = this%nodes           
           do e1 = 3, this%degMax
              this%listAdj(this%aux(e1)) = 1           
              this%listAdj(this%aux(e1) + 1) = e1 - 1
              this%listAdj(this%aux(e1) + 2) = e1 + 1              
           enddo
           this%listAdj(this%aux(this%nodes)) = 1
           this%listAdj(this%aux(this%nodes) + 1) = this%degMax
           this%listAdj(this%aux(this%nodes) + 2) = 2
           !------------------------------------------------------------
           ! /Para os raios da roda
           !------------------------------------------------------------                      
        end subroutine   
!-----------------------------------------------------------------------
        subroutine liga_grafo_roda_sob_ataque(this, grau_hub, p_at, sem10)
           class(grafo_roda) :: this
           type(rndgen) :: gen10
           integer, intent(in) :: grau_hub
           integer :: e1, e12, e2, e23, e3
           integer :: cont
           real(dp), intent(in) :: p_at
           integer :: sem10
           real(dp) :: prob
           integer, allocatable :: listAdjAux(:), degAux(:), AuxAux(:)
           !------------------------------------------------------------           
           call gen10%init(sem10)
           !------------------------------------------------------------           
           call this%liga_roda(grau_hub)
           !if1---------------------------------------------------------
           if( p_at > 0)then
              !---------------------------------------------------------
              if(allocated(degAux)) deallocate(degAux)
              allocate(degAux(size(this%deg)))
              !---------------------------------------------------------
              do e1 = 1, this%nodes
                 degAux(e1) = this%deg(e1)
              enddo
              !---------------------------------------------------------           
              do e1 = 2, this%nodes
                 do e12 = this%aux(e1), this%aux(e1) + this%deg(e1) - 1
                    e2 = this%listAdj(e12)
                    if( e2 == 1 ) cycle
                    if(this%listAdj(e12) == 0) cycle
                    prob = gen10%rnd()
                    !---------------------------------------------------
                    if(prob <= p_at)then                       
                       this%listAdj(e12) = 0
                       !------------------------------------------------
                       degAux(e1) = degAux(e1) - 1
                       !------------------------------------------------
                       do e23 = this%aux(e2), this%aux(e2) + this%deg(e2) - 1
                          e3 = this%listAdj(e23)
                          if( e3 /= e1) cycle
                          this%listAdj(e23) = 0
                          !---------------------------------------------
                          degAux(e2) = degAux(e2) - 1
                          !---------------------------------------------
                       enddo
                       !------------------------------------------------
                    endif
                    !---------------------------------------------------                    
                 enddo
              enddo
              !---------------------------------------------------------
              if(allocated(AuxAux)) deallocate(AuxAux)
              allocate(AuxAux(size(this%aux)))
              AuxAux(1) = 1
              do e1 = 2, this%nodes
                 AuxAux(e1) = AuxAux(e1-1) + degAux(e1-1)
              enddo
              !---------------------------------------------------------
              if(allocated(listAdjAux)) deallocate(listAdjAux)
              allocate(listAdjAux(sum(degAux)))
              !---------------------------------------------------------
              degAux = 0
              !---------------------------------------------------------           
              do e1 = 1, this%nodes
                 do e12 = this%aux(e1), this%aux(e1) + this%deg(e1) - 1
                    e2 = this%listAdj(e12)
                    if(e2 == 0) cycle
                    listAdjAux(AuxAux(e1) + degAux(e1)) = e2
                    degAux(e1) = degAux(e1) + 1 
                 enddo
              enddo
              !---------------------------------------------------------
              deallocate(this%listAdj)
              allocate(this%listAdj(size(listAdjAux)))
              !---------------------------------------------------------
              do e1 = 1, size(listAdjAux)
                 this%listAdj(e1) = listAdjAux(e1)
              enddo
              !---------------------------------------------------------
              deallocate(listAdjAux)
              !---------------------------------------------------------
              do e1 = 1, size(this%aux)
                 this%aux(e1) = AuxAux(e1)
              enddo
              !---------------------------------------------------------
              deallocate(AuxAux)
              !---------------------------------------------------------
              do e1 = 1, size(this%deg)
                 this%deg(e1) = degAux(e1)
              enddo
              !---------------------------------------------------------
              deallocate(degAux)
           !------------------------------------------------------------              
           endif
           this%degMin = minval(this%deg)
           this%degMax = maxval(this%deg)
           this%degMean = sum(this%deg)/(1.0d0 * this%nodes)
           !------------------------------------------------------------
           !endif1------------------------------------------------------
        end subroutine
!-----------------------------------------------------------------------
        subroutine liga_grafo_estrela(this, grau_estrela)
           class(grafo_estrela) :: this
           integer, intent(in) :: grau_estrela
           integer :: e1, e2
           
           this%nodes = grau_estrela + 1
           
           if(allocated(this%deg)) deallocate(this%deg)
              allocate(this%deg(this%nodes))
           
           this%deg(1) = grau_estrela
           
           do e1 = 2, this%nodes
              this%deg(e1) = 1   
           enddo   
           this%degMin = 1
           this%degMax = grau_estrela
           this%sumDeg = 2 * grau_estrela
           this%degMean = 1.0_dp * this%sumDeg/(1.0_dp * this%nodes)
           
           if(allocated(this%aux)) deallocate(this%aux)
              allocate(this%aux(this%nodes))           
           
           this%aux(1) = 1
           
           do e1 = 2, this%nodes
              this%aux(e1) = this%aux(e1-1) + this%deg(e1-1)
           enddo

           if(allocated(this%listAdj)) deallocate(this%listAdj)
              allocate(this%listAdj(this%sumDeg))           
           
           do e1 = 1, this%degMax
              this%listAdj(e1) = e1 + 1
              this%listAdj(this%degMax + e1) = 1
           enddo
           
        end subroutine   

		subroutine ligaRedeReal(this,arquivo, label)
			class(grafo) :: this
			character(len=*) :: arquivo
			integer, intent(in) :: label
			integer :: node1, node2, arestas, N_prov
			integer :: error
			integer :: i
			
			arquivo = trim(adjustl(arquivo))
			
			open(unit=label, file=arquivo, status='old')
			
			arestas = 0
			N_prov = 0
			
!#######################################################################			
			do
				read(label,*,iostat=error) node1, node2
				if(error /= 0) exit
				arestas = arestas + 1
				N_prov = max(node1, N_prov)
				N_prov = max(node2, N_prov)
			enddo
			rewind(label)
			write(*,*) 'numero de arestas', arestas
			write(*,*) "Maior indice do sitio da rede eh ", N_prov
!#######################################################################			
			
			call this%iniciaGrafo(N_prov)
			
			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(arestas,2))
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(2*arestas))			
			
			this%edge = arestas
			
!#######################################################################			
			do i = 1, arestas
				read(label,*) node1, node2
				this%deg(node1) = this%deg(node1) + 1
				this%deg(node2) = this%deg(node2) + 1
				this%matriz(i, 1) = node1
				this%matriz(i, 2) = node2
			enddo			
!##########################################################						
			
			this%aux(1) = 1
			
			do i = 2, this%nodes
				this%aux(i) = this%aux(i-1) + this%deg(i-1)
			enddo

!##########################################################
			
			this%deg = 0
			
			do i = 1, this%edge
				this%listAdj(this%aux(this%matriz(i, 1)) + this%deg(this%matriz(i,1))) = this%matriz(i,2)
					this%deg(this%matriz(i,1)) = this%deg(this%matriz(i,1)) + 1
				this%listAdj(this%aux(this%matriz(i, 2)) + this%deg(this%matriz(i,2))) = this%matriz(i,1)				
					this%deg(this%matriz(i,2)) = this%deg(this%matriz(i,2)) + 1
			enddo
	
			this%degMin = minval(this%deg)
			this%degMax = maxval(this%deg)
			
			this%degMean = 1.0_dp * sum(this%deg)/this%nodes
			
		end subroutine
		
!#######################################################################		
		!###############################################################
		!	Grafo Scale-free rico em folhas
		!###############################################################
		
		subroutine inicia_PL_UCM_Folhas(this, prob_folhas2, k_Min2, k_Max2, gama2, semente)
			class(grafo_PL_UCM_Folhas) :: this
			integer, intent(in) :: k_Min2, k_Max2
			real(dp), intent(in) :: prob_folhas2, gama2
			integer :: semente

			real(dp), allocatable :: D_k(:)
		
			real(dp) :: prob, sum_Dk, A1
			integer :: i1, i2, i3, k1, the_chosen
			type(rndgen) :: gen
			
			!###########################################################
			!	Inicializacao de variaveis importantes e atributos 
			!	de classe
			!###########################################################
			
			this%gamm = gama2

			this%prob_folhas = prob_folhas2			
			!###########################################################
			!	Alocamento da lista auxiliar
			!###########################################################
			
			if(allocated(D_k)) deallocate(D_k)
				allocate(D_k(k_Min2:k_Max2))


			!###########################################################
			!	Inicializacao da lista auxiliar
			!###########################################################			

				D_k(k_Max2) = (1.0_dp * k_Max2)**(-this%gamm)				
				
				do k1 = k_Max2 - 1, k_Min2, -1
					D_k(k1) = (1.0_dp * k1) **(-this%gamm)
					D_k(k1) = D_k(k1) + D_k(k1 + 1)
				enddo

			!###########################################################
			!	Normalizacao da lista.
			!###########################################################	
					sum_Dk = D_k(k_Min2)
			!###########################################################
			!	Normalizacao da lista auxiliar
			!	e calculo da distribuicao com tubos.
			!###########################################################

				D_k = D_k/sum_Dk
				

			!###########################################################
			!	Criacao da lista de graus
			!###########################################################				
				
				call gen%init(semente)
				
loop_noh1:		do i1 = 1, this%nodes
					prob = gen%rnd()
					
					if(prob > this%prob_folhas)then
						prob = gen%rnd()
	loop_grau1:			do i2 = k_Max2, k_Min2, -1
							
							if(prob <= D_k(i2)) then
									this%deg(i1) = i2
									exit loop_grau1
							endif
						enddo	loop_grau1
					else
						this%deg(i1) = 1
					endif	
				enddo loop_noh1
				
			!###########################################################
			!	Redefinimos degMax de acordo com a distribuicao
			!	obtida.			
			!###########################################################
			
					
!#######################################################################
!	Se o numero total de stubs nao eh par, eu somo 1 unidade a um
!	no arbitrario da rede.
!#######################################################################
			this%sumDeg = sum(this%deg)
						
			if(mod(this%sumDeg, 2)/=0) then
				the_chosen = this%nodes
				this%deg(the_chosen) = this%deg(the_chosen) + 1
				this%sumDeg = this%sumDeg + 1
			endif		
			deallocate(D_k)


!#######################################################################
!	Refaco essa parte, porque o grau min e o grau maximo podem mudar.
!	Se eu disser que o grau_min > 2, no final das contas,
!	o grau_min = 2. Se ele for 1, o grau min e 1 sim.
!	No entanto, o grau maximo pode ser determinado pela estatistica
!	e ser menor do que o estipulado.
!#######################################################################			
			this%degMin = minval(this%deg)
			this%degMax = maxval(this%deg)
!#######################################################################
!	Finalmente aloca as ultimas listas atributos de classe
!	do objeto.
!#######################################################################
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
			
				
			this%aux(1) = 1
			
			do i1 = 2, this%nodes
				this%aux(i1) = this%aux(i1-1) + this%deg(i1-1)
			enddo
							
		end subroutine
	
	
!#######################################################################
		
		subroutine iniciaGrafo(this, N)
			integer, intent(in) :: N
			class(grafo), intent(inout) :: this
		
			this%nodes = N
		
			if(allocated(this%deg)) deallocate(this%deg)
				allocate(this%deg(N))
				this%deg = 0			
			if(allocated(this%aux)) deallocate(this%aux)
				allocate(this%aux(N))				
			
			!###########################################################
			! Nao tem como iniciar matriz, aux e nem listAdj
			! ateh que se conheca a distribuicao de graus
			!###########################################################
			
			
		end subroutine
		
		
		
	
		
!#######################################################################
!	Modelo UCM para Rede Scale Free
!#######################################################################        
		subroutine inicia_PL_UCM(this, k_Min2, k_Max2, gama2, semente)

!#######################################################################
!	Argumentos da subrotina
!#######################################################################

			class(grafo_PL_UCM) :: this
			integer, intent(in) :: k_Min2
			real(dp), intent(in) :: k_Max2
			real(dp), intent(in) :: gama2
			integer :: semente
					
!#######################################################################
!	Lista auxiliar	
!#######################################################################			

			real(dp), allocatable :: D_k(:)

!#######################################################################
!	Variaveis auxiliares
!#######################################################################

			real(dp) :: prob, sum_Dk, harvest
			integer :: i1, i2, i3, k1, the_chosen
			integer :: tic, toc, taxa
			type(rndgen) :: gen
					
			
			
			!###########################################################
			!	Inicializacao de variaveis importantes e atributos 
			!	de classe
			!###########################################################
			this%gamm = gama2
			this%degMin = k_Min2
			this%degMax = k_Max2
			
			!###########################################################
			!	Alocamento da lista auxiliar
			!###########################################################
			
			if(allocated(D_k)) deallocate(D_k)
				allocate(D_k(this%degMin:this%degMax))


			!###########################################################
			!	Inicializacao da lista auxiliar
			!###########################################################			

				D_k(this%degMax) = (1.0_dp * this%degMax)**(-this%gamm)				
				
				sum_Dk = D_k(this%degMax)
				
				do k1 = this%degMax - 1, this%degMin, -1
					D_k(k1) = (1.0_dp * k1) **(-this%gamm)
					D_k(k1) = D_k(k1) + D_k(k1 + 1)
				enddo

				sum_Dk = D_k(this%degMin)

			!###########################################################
			!	Normalizacao da lista auxiliar
			!###########################################################

				D_k = D_k/sum_Dk
				

			!###########################################################
			!	Criacao da lista de graus
			!###########################################################				
				call gen%init(semente)
				
loop_noh:		do i1 = 1, this%nodes
					prob = gen%rnd()
					
loop_grau:         do i2 = this%degMin, this%degMax-1
                      if(prob > D_k(i2+1))then
                         this%deg(i1) = i2
                         exit loop_grau
                      elseif(prob <= D_k(this%degMax))then
                          this%deg(i1) = this%degMax
                          exit loop_grau
                      endif
                   enddo loop_grau
                   
					
				enddo loop_noh
!#######################################################################
!	Se o numero total de stubs nao eh par, eu somo 1 unidade a um
!	no tomado aleatoriamente da rede. Como os hubs devem ser 
!	poucos comparado com os demais nos, nao ha preocupacao em que
!	um hub seja escolhido
!#######################################################################
			this%sumDeg = sum(this%deg)
						

			if(mod(this%sumDeg, 2)/=0) then
				the_chosen = this%nodes
				this%deg(the_chosen) = this%deg(the_chosen) + 1
				this%sumDeg = this%sumDeg + 1
			endif		
			deallocate(D_k)

			this%degMax = maxval(this%deg)
!#######################################################################
!	Finalmente aloca as ultimas listas atributos de classe
!	do objeto.
!#######################################################################
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
			allocate(this%listAdj(this%sumDeg))
							
			this%aux(1) = 1
			
			do i1 = 2, this%nodes
				this%aux(i1) = this%aux(i1-1) + this%deg(i1-1)
			enddo
				
		end subroutine




!#######################################################################
!	Rede gerada por modelo UCM com tubos
!#######################################################################


		subroutine inicia_PL_UCM_Tubes(this, prob_tube2, k_Min2, k_Max2, gama2, semente)

!#######################################################################
!	Argumentos da subrotina
!#######################################################################

			class(grafo_PL_UCM_Tubes) :: this
			integer, intent(in) :: k_Min2, k_Max2
			real(dp), intent(in) :: prob_tube2, gama2
			integer :: semente
					
!#######################################################################
!	Lista auxiliar	
!#######################################################################			

			real(dp), allocatable :: D_k(:)

!#######################################################################
!	Variaveis auxiliares
!#######################################################################

			real(dp) :: prob, sum_Dk, A1
			integer :: i1, i2, i3, k1, the_chosen
			type(rndgen) :: gen
					
			

			!###########################################################
			!	Inicializacao de variaveis importantes e atributos 
			!	de classe
			!###########################################################
			this%gamm = gama2
			this%degMin = 2
			this%degMax = k_Max2
			this%prob_tube = prob_tube2			
			!###########################################################
			!	Alocamento da lista auxiliar
			!###########################################################
			
			if(allocated(D_k)) deallocate(D_k)
				allocate(D_k(this%degMin:this%degMax))


			!###########################################################
			!	Inicializacao da lista auxiliar
			!###########################################################			

				D_k(this%degMax) = (1.0_dp * this%degMax)**(-this%gamm)				
				
				do k1 = this%degMax - 1, this%degMin + 1, -1
					D_k(k1) = (1.0_dp * k1) **(-this%gamm)
					D_k(k1) = D_k(k1) + D_k(k1 + 1)
				enddo

			!###########################################################
			!	Normalizacao da lista.
			!###########################################################	
					sum_Dk = D_k(this%degMin + 1)
			!###########################################################
			!	Normalizacao da lista auxiliar
			!	e calculo da distribuicao com tubos.
			!###########################################################

				D_k = D_k/sum_Dk
				

			!###########################################################
			!	Criacao da lista de graus
			!###########################################################				
				
				call gen%init(semente)
				
loop_noh1:		do i1 = 1, this%nodes
					prob = gen%rnd()
					
					if(prob > this%prob_tube)then
						prob = gen%rnd()
	loop_grau1:			do i2 = this%degMax, this%degMin+1, -1
							
							if(prob <= D_k(i2)) then
									this%deg(i1) = i2
									exit loop_grau1
							endif
						enddo	loop_grau1
					else
						this%deg(i1) = 2
					endif	
				enddo loop_noh1
				
			!###########################################################
			!	Redefinimos degMax de acordo com a distribuicao
			!	obtida.			
			!###########################################################
			
					
!#######################################################################
!	Se o numero total de stubs nao eh par, eu somo 1 unidade a um
!	no arbitrario da rede.
!#######################################################################
			this%sumDeg = sum(this%deg)
						
			if(mod(this%sumDeg, 2)/=0) then
				the_chosen = this%nodes
				this%deg(the_chosen) = this%deg(the_chosen) + 1
				this%sumDeg = this%sumDeg + 1
			endif		
			deallocate(D_k)


!#######################################################################
!	Refaco essa parte, porque o grau min e o grau maximo podem mudar.
!	Se eu disser que o grau_min > 2, no final das contas,
!	o grau_min = 2. Se ele for 1, o grau min e 1 sim.
!	No entanto, o grau maximo pode ser determinado pela estatistica
!	e ser menor do que o estipulado.
!#######################################################################			
			this%degMin = minval(this%deg)
			this%degMax = maxval(this%deg)
!#######################################################################
!	Finalmente aloca as ultimas listas atributos de classe
!	do objeto.
!#######################################################################
			
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
			
				
			this%aux(1) = 1
			
			do i1 = 2, this%nodes
				this%aux(i1) = this%aux(i1-1) + this%deg(i1-1)
			enddo
				
		end subroutine


		!###############################################################
		!		Subrotina que conecta os nos do grafo RRN
		!###############################################################
	
	
		subroutine ligaUCM(this, seed, escreveListaConexoes)
!#######################################################################
!		objeto tipo grafo
!#######################################################################
			class(grafo) :: this

!#######################################################################
!		Lista auxiliar contendo 'stubs'
!#######################################################################
			integer, allocatable :: lC(:)

!#######################################################################
!		Lista auxiliar genericas
!#######################################################################
			integer :: i, j, k
			integer :: seed, seed2
			integer :: cand1, cand2, N1
			integer :: lastC, nStubs, stub1, stub2, sumDeg2
			integer (kind=8) :: conta_falhas1, &		    
			conta_falhas2, teto_de_falhas1, teto_de_falhas2
			
!#######################################################################
!		Lista auxiliar para ir contando o grau, que cresce.
!		Tudo indica que esta lista sera deprecada no codigo.
!#######################################################################			
			integer, allocatable :: degAux(:), lista_auxiliar_limpeza(:)

!#######################################################################
!		Variavel logica que vai dizer se auto-ligacoes
!		e ligacoes multiplas foram evitadas
!#######################################################################						
			logical :: ligouStubs
			
!#######################################################################
!		Objeto do tipo gerador de numeros pseudo aleatorios
!#######################################################################			
			type(rndgen) :: gen
			
			real(dp) :: esquenta_gerador

!#######################################################################
!	Decide se escreve arquivo com lista de conexoes
!#######################################################################			
			logical :: escreveListaConexoes

!#######################################################################
!		Inicializamos o gerador de numeros pseudo-aleatorios
!		e em seguida o esquentamos.
!#######################################################################
			call gen%init(seed)

			this%edge=0
			
			this%sumDeg = sum(this%deg)

!#######################################################################
!	Se gera_lista_Adjacencia = .True., aloca a matriz
!#######################################################################
			

			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%sumDeg/2, 2))


!#######################################################################
!	Gera lista auxiliar.
!#######################################################################
			
			if(allocated(lC)) deallocate(lC)
				allocate(lC(this%sumDeg))


!#######################################################################
!	Preenche a lista auxiliar
!#######################################################################
							
			k = 1
			do i = 1, this%nodes
				if(this%deg(i) == 0) cycle
				do j = 1, this%deg(i)
					lC(k) = i
					k = k + 1
				enddo
			enddo

!#######################################################################
!	Checa os stubs disponiveis
!#######################################################################
			
			lastC = size(lC)
			nStubs = lastC/2
!			write(*,*) "numero de stubs disponiveis", nStubs


!#######################################################################
!	Outra lista auxiliar
!#######################################################################
			
			if(allocated(degAux)) deallocate(degAux)
			allocate(degAux(this%nodes))
			degAux = 0
			
			
!#######################################################################
!		Antes o loop era orientado aos nos, agora
!		o farei orientado a stubs, como parece ser sugerido
!		na literatura.
!#######################################################################			

conta_falhas2=0		
loop_stubs: do while (lastC > 0)
!		
		stub1 = gen%int(1, lastC)
		cand1 = lC(stub1)

!#######################################################################
!		Loop que testa se auto-conexoes foram evitadas.
!		LigaStubs = .False., porque somos pessimistas,
!		porem pagamos pra ver qual eh.
!#######################################################################			
		ligouStubs = .False.

!#######################################################################
!	Criei um contador de falhas de conexao para um stub1 arbitrario.
!#######################################################################
		
		conta_falhas1 = 0

!#######################################################################
!	Criei um contador de falhas para stubs em geral.
!	Se um teto de falhas for alcancado, interrompe o processo.
!#######################################################################	
		
		teto_de_falhas2 = 10 * lastC
		if(conta_falhas2 > teto_de_falhas2)then
!			write(*,*) " "
!			write(*,*) " Teremos que interromper a rotina com ", conta_falhas2, " falhas no primeiro stub."
!			write(*,*) " "
!			write(*,*) "Temos ", lastC, " stubs disponiveis ainda e tivemos que interromper"
!			write(*,*) " "
!			write(*,*) "Foram feitas ", this%edge, " ligacoes."
!			write(*,*) " "
!			write(*,*) " Interrompendo ligacoes e passando para a reestruturacao das listas..."
!			write(*,*) " "
			exit loop_stubs
		endif
Principal:		do while(.not. ligouStubs)
						
						
						!###############################################
						!	Se se permanece nesse loop por um
						!	tempo maior que o maior que
						!	100 * o numero de stubs disponiveis,
						!	terminamos o processo
						!	Pus uma condicao de parada meio ousada
						!	aqui em baixo. Em observacao!!!
						!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!###############################################

						teto_de_falhas1 = 10 * lastC
						if(conta_falhas1 > teto_de_falhas1)then
							conta_falhas2 = conta_falhas2 + 1
							exit Principal
						endif
						
						stub2 = gen%int(1, lastC)
						cand2 = lC(stub2)
						
!#######################################################################
!		Toda vez que um no se liga, seu degAux aumenta em uma unidade
!		e um de seus stubs some de lC.
!		Quando seu degAux == this%deg, ele ja nao estara mais em lC.
!		Portanto, essa condicao abaixo nunca sera satisfeita.

!						if(degAux(cand2) == this%deg(cand2))then
!							conta_falhas1 = conta_falhas1 + 1
!							cycle Principal
!						endif
!
!		e nem essa
!
!						if(degAux(cand1) == this%deg(cand1)) &
!						exit Principal
!#######################################################################





						
!#######################################################################
!		Aqui evito autoconexoes e ligacoes repetidas.
!		Se nao ha ligacoes repetidas apos o 'else', ligaStubs
!		muda de valor, permitindo sair do loop .
!		Toda vez que a conexao falha, contamos para entao resetar
!		a semente, caso ocorra um numero muito grande de falhas.
!#######################################################################
						
						if(this%deg(cand1) == 1)then
							if(this%deg(cand2) == 1)then
								conta_falhas1 = conta_falhas1 + 1
								cycle Principal
							endif
						endif
						
						if(cand2 == cand1)then
							conta_falhas1 = conta_falhas1 + 1
							cycle Principal
						else
							do k = this%aux(cand1), this%aux(cand1) + degAux(cand1) - 1
								if(this%listAdj(k) == cand2) then
									conta_falhas1 = conta_falhas1 + 1
									cycle Principal
								endif
							enddo
							ligouStubs = .True.
						endif
										
						this%listAdj(this%aux(cand1)+degAux(cand1)) = cand2
						this%listAdj(this%aux(cand2)+degAux(cand2)) = cand1
						degAux(cand1) = degAux(cand1) + 1
						degAux(cand2) = degAux(cand2) + 1
						
						this%edge = this%edge + 1

						conta_falhas1 = 0

!#######################################################################
!		Se escolhi gerar a lista de adjacencias pro Gephi,
!		ele gera, caso contrario, ignora.
!#######################################################################
						

						this%matriz(this%edge, 1) = cand1
						this%matriz(this%edge, 2) = cand2

!#######################################################################
!	Se stub2 vier depois de stub1, nao podemos eliminar o stub1
!	primeiro pois, se stub2 for lastC, colocaremos cand2 (stub2)
!	na posicao stub1 e, quando formos eliminar o stub2, ele estara
!	virtualmente fora da lista lC. Dito isso, ao colocarmos o
!	lastC atual na posicao stub2, estaremos tirando um provavel lastC,
!	que nao foi ligado ainda, da lista lC.
!	Situacao similar ocorre se pensarmos inicialmente no stub1.
!#######################################################################

                        if(stub1 == lastC)then
                           lastC = lastC - 1                           
                           lC(stub2) = lC(lastC)
                           lastC = lastC - 1
                        elseif(stub2 == lastC)then
                           lastC = lastC - 1                           
                           lC(stub1) = lC(lastC)
                           lastC = lastC - 1
                        else 
							lC(stub1) = lC(lastC)
							lastC = lastC - 1
							lC(stub2) = lC(lastC)
							lastC = lastC - 1
						endif
			enddo Principal
	enddo loop_stubs
!#######################################################################
!	Devido a possibilidade de o modelo UCM poder deixar stubs soltos,
!	vamos 'confirmar' a rede. Neste caso, criar a matriz de conexoes
!	passa a nao ser mais opcional, mas obrigatorio.
!
!#######################################################################			


			
!			write(*,*) "Sobraram ", nStubs - this%edge," stubs"

!#######################################################################
!	Se sobraram stubs, vamos reformular a lista de graus.
!	A lista e zerada e a partir da matriz de adjacencias,
!	nos recontamos os graus.
!#######################################################################
			if((nStubs - this%edge) > 0)then

			!###########################################################
			!	Se o numero de conexoes nao eh o previsto,
			!	zeramos a lista de graus e recontamos a lista de
			!	conexoes atraves da matriz de conexoes.
			!###########################################################
				
			   write(*,*) " "
			   write(*,*) "Faremos a reestruturacao das listas da rede"
			   write(*,*) " "
			   write(*,*) "Parametros: "
			   write(*,*) " "
			   write(*,*) "Numero de nos: ", this%nodes
			   write(*,*) "Numero de ligacoes: ", this%edge
			   write(*,*) "Grau maximo: ", this%degMax
			   write(*,*) "Grau medio: ", this%degMean 
			   write(*,*) " "
		

               this%deg = 0
				
               do i = 1, this%edge
                  cand1 = this%matriz(i,1)
                  this%deg(cand1) = this%deg(cand1) + 1
					
                  cand2 =	this%matriz(i,2)
                  this%deg(cand2) = this%deg(cand2) + 1
               enddo

			!###########################################################
			!	A soma de graus e recontada
			!###########################################################			
				

!				write(*,*) " "
!				write(*,*) "Contamos ", this%sumDeg, " stubs usados em ligacoes na rede ..."
!				write(*,*) " "

			!###########################################################
			!	Guardamos o numero original de nos para testar se houve
			!	discrepancia.
			!###########################################################			
				
				N1 = this%nodes
				

			!###########################################################
			!	Zeramos o numero de nos original e recontamos a partir da lista
			!	de graus
			!###########################################################			
				this%nodes = 0
				do i = 1, N1
					if(this%deg(i) > 0) this%nodes = this%nodes + 1
				enddo

!				write(*,*) " "
!				write(*,*) "Contamos ", this%nodes, " nos ligantes ate agora..."
!				write(*,*) " "

!#######################################################################
!	Se o numero de nos decaiu, desalocamos e alocamos a lista
!	de graus e a lista auxiliar e recontamos a lista de graus.
!#######################################################################			

! Minha this%matriz ainda guarda os indices originais dos nos.
! Quando this%nodes < N1, da merda.
! Minha ideia eh reindexar os nos do de indice menor para o de indice menor.
! Mas como?


				
				if(this%nodes < N1)then

!#######################################################################			
!	Mapeamento dos indices antigos para os atuais

					if(allocated(lista_auxiliar_limpeza)) deallocate(lista_auxiliar_limpeza)
						allocate(lista_auxiliar_limpeza(N1))
						lista_auxiliar_limpeza = 0
					
					do i = 1, this%edge
						lista_auxiliar_limpeza(this%matriz(i, 1)) = 1
						lista_auxiliar_limpeza(this%matriz(i, 2)) = 1
					enddo	
					
					do i = N1, 1, -1
						if(lista_auxiliar_limpeza(i) > 0 )then
							lista_auxiliar_limpeza(i) = this%nodes
							this%nodes = this%nodes - 1
						endif  
					enddo	
					
					do i = 1, this%edge
						cand1 = this%matriz(i,1)
						this%matriz(i,1) = lista_auxiliar_limpeza(cand1)
						cand2 = this%matriz(i,2)
						this%matriz(i,2) = lista_auxiliar_limpeza(cand2)
					enddo
					
					this%nodes = maxval(this%matriz)
!#######################################################################						
							
					if(allocated(this%deg)) deallocate(this%deg)
					allocate(this%deg(this%nodes))
					this%deg = 0
						
					if(allocated(this%aux)) deallocate(this%aux)
					allocate(this%aux(this%nodes))
					
					do i = 1, this%edge
						cand1 =	this%matriz(i,1)
						this%deg(cand1) = this%deg(cand1) + 1

						cand2 =	this%matriz(i,2)
						this%deg(cand2) = this%deg(cand2) + 1
					enddo
					
				endif

!#######################################################################
!	Recontamos a lista auxiliar.
!#######################################################################			

				
				this%aux(1) = 1
					
				do i = 2, this%nodes
					this%aux(i) = this%aux(i-1) + this%deg(i-1)
				enddo

!#######################################################################
!	Finalmente, recontamos a lista de adjacencias.
!#######################################################################			

				this%sumDeg = sum(this%deg)

				if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
				this%deg = 0
				
!#######################################################################
!	Para isto, re-zeramos a lista de graus e comecamos a contar
!	tudo de novo... isso ate achar um metodo melhor.
!#######################################################################							
				do i = 1, this%edge
					
					cand1 =	this%matriz(i,1)
					cand2 =	this%matriz(i,2)
					this%listAdj(this%aux(cand1) + this%deg(cand1)) = cand2
					
					this%deg(cand1) = this%deg(cand1) + 1
					
					this%listAdj(this%aux(cand2) + this%deg(cand2)) = cand1
					
					this%deg(cand2) = this%deg(cand2) + 1					
				enddo
			endif

			this%degMax = maxval(this%deg)
			this%degMin = minval(this%deg)
			
			write(*,*) 'O grau maximo apos a ligacao da rede eh : ', this%degMax
			write(*,*) 'O grau minimo apos a ligacao da rede eh : ', this%degMin
						
			if(escreveListaConexoes)then
				open(113,file='lista_conexoes.csv', status='unknown')

				do i = 1, this%edge
					write(113,*) this%matriz(i,1), ",", this%matriz(i,2)
				enddo
				close(113)
			endif
			
			this%degMean = 1.0_dp * this%sumDeg/this%nodes
			
			write(*,*) 'Este e o grau medio apos a ligacao da rede: ', this%degMean
	end subroutine		




!#######################################################################

	
		subroutine religaUCM_n_imune(this, seed, escreveListaConexoes)
!#######################################################################
!		objeto tipo grafo
!#######################################################################
			class(grafo) :: this

!#######################################################################
!		Lista auxiliar contendo 'stubs'
!#######################################################################
			integer, allocatable :: lC(:), lista_auxiliar_limpeza(:)

!#######################################################################
!		Lista auxiliar genericas
!#######################################################################
			integer :: i, j, k
			integer :: seed, seed2
			integer :: cand1, cand2, N1
			integer :: lastC, nStubs, stub1, stub2, sumDeg2
			integer (kind=8) :: conta_falhas1,&
			conta_falhas2, teto_de_falhas1, teto_de_falhas2
			
!#######################################################################
!		Lista auxiliar para ir contando o grau, que cresce.
!		Tudo indica que esta lista sera deprecada no codigo.
!#######################################################################			
			integer, allocatable :: degAux(:)

!#######################################################################
!		Variavel logica que vai dizer se auto-ligacoes
!		e ligacoes multiplas foram evitadas
!#######################################################################						
			logical :: ligouStubs
			
!#######################################################################
!		Objeto do tipo gerador de numeros pseudo aleatorios
!#######################################################################			
			type(rndgen) :: gen
			
			real(dp) :: esquenta_gerador

!#######################################################################
!	Decide se escreve arquivo com lista de conexoes
!#######################################################################			
			logical :: escreveListaConexoes

!#######################################################################
!		Inicializamos o gerador de numeros pseudo-aleatorios
!		e em seguida o esquentamos.
!#######################################################################
			integer :: conta_sobras


!#######################################################################
! \	Aqui vai fazer toda a diferenca
!	O resto fora daqui eh identico
!#######################################################################
			conta_sobras = 0

!#######################################################################
!	Nesta listra grau2, so entram nos que fazem parte da componente
!	gigante.
!#######################################################################
				
			do i =1, this%nodes
				if(grau2(i) > 0) conta_sobras = conta_sobras + 1
			enddo
			
			if(allocated(this%deg)) deallocate(this%deg)
				allocate(this%deg(conta_sobras))
				this%deg = 0
			
			k = 0
			do i =1, this%nodes
				if(grau2(i) > 0)then
					k = k + 1
					this%deg(k) = grau2(i)
				endif	
			enddo
			
			this%nodes = conta_sobras
							
				
!#######################################################################
!	Atualiza propriedades associadas a this%deg
!#######################################################################
			
			this%sumDeg = sum(this%deg)

			this%degMax = maxval(this%deg)

!#######################################################################
!
!#######################################################################

			if(mod(this%sumDeg,2) > 0)then
				i = gen%int(1, this%nodes)
				this%deg(i) = this%deg(i) + 1
				
				this%degMax = max(this%degMax, this%deg(i))
				
				this%sumDeg = this%sumDeg + 1
				
			endif
			
			
			
			if(allocated(this%aux)) deallocate(this%aux)
				allocate(this%aux(this%nodes))
			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
							
				this%aux(1) = 1
				do i = 2, this%nodes
					this%aux(i) = this%aux(i-1) + this%deg(i-1)
				enddo		
!#######################################################################
!	Aqui vai fazer toda a diferenca /
!#######################################################################
			
	
			call gen%init(seed)
			
			do i = 1, 10 * this%degMax
				esquenta_gerador = gen%rnd()
			enddo
			
			this%edge=0
			

!#######################################################################
!	Se gera_lista_Adjacencia = .True., aloca a matriz
!#######################################################################
			

			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%sumDeg/2, 2))


!#######################################################################
!	Gera lista auxiliar.
!#######################################################################
			
			if(allocated(lC)) deallocate(lC)
				allocate(lC(this%sumDeg))


!#######################################################################
!	Preenche a lista auxiliar
!#######################################################################
							
			k = 1
			do i = 1, this%nodes
				if(this%deg(i) == 0) cycle
				do j = 1, this%deg(i)
					lC(k) = i
					k = k + 1
				enddo
			enddo

!#######################################################################
!	Checa os stubs disponiveis
!#######################################################################
			
			lastC = size(lC)
			nStubs = lastC/2
!			write(*,*) "numero de stubs disponiveis", nStubs


!#######################################################################
!	Outra lista auxiliar
!#######################################################################
			
			if(allocated(degAux)) deallocate(degAux)
			allocate(degAux(this%nodes))
			degAux = 0
			
			
!#######################################################################
!		Antes o loop era orientado aos nos, agora
!		o farei orientado a stubs, como parece ser sugerido
!		na literatura.
!#######################################################################			

conta_falhas2=0		
loop_stubs: do while (lastC > 0)
!		
		stub1 = gen%int(1, lastC)
		cand1 = lC(stub1)

!#######################################################################
!		Loop que testa se auto-conexoes foram evitadas.
!		LigaStubs = .False., porque somos pessimistas,
!		porem pagamos pra ver qual eh.
!#######################################################################			
		ligouStubs = .False.

!#######################################################################
!	Criei um contador de falhas de conexao para um stub1 arbitrario.
!#######################################################################
		
		conta_falhas1 = 0

!#######################################################################
!	Criei um contador de falhas para stubs em geral.
!	Se um teto de falhas for alcancado, interrompe o processo.
!#######################################################################	
		
		teto_de_falhas2 = 10 * lastC
		if(conta_falhas2 > teto_de_falhas2)then
!			write(*,*) " "
!			write(*,*) " Teremos que interromper a rotina com ", conta_falhas2, " falhas no primeiro stub."
!			write(*,*) " "
!			write(*,*) "Temos ", lastC, " stubs disponiveis ainda e tivemos que interromper"
!			write(*,*) " "
!			write(*,*) "Foram feitas ", this%edge, " ligacoes."
!			write(*,*) " "
!			write(*,*) " Interrompendo ligacoes e passando para a reestruturacao das listas..."
!			write(*,*) " "
			exit loop_stubs
		endif
Principal:		do while(.not. ligouStubs)
						
						
						!###############################################
						!	Se se permanece nesse loop por um
						!	tempo maior que o maior que
						!	100 * o numero de stubs disponiveis,
						!	terminamos o processo
						!	Pus uma condicao de parada meio ousada
						!	aqui em baixo. Em observacao!!!
						!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!###############################################

						teto_de_falhas1 = 10 * lastC
						if(conta_falhas1 > teto_de_falhas1)then
							conta_falhas2 = conta_falhas2 + 1
							exit Principal
						endif
						
						stub2 = gen%int(1, lastC)
						cand2 = lC(stub2)
						
!#######################################################################
!		Toda vez que um no se liga, seu degAux aumenta em uma unidade
!		e um de seus stubs some de lC.
!		Quando seu degAux == this%deg, ele ja nao estara mais em lC.
!		Portanto, essa condicao abaixo nunca sera satisfeita.

!						if(degAux(cand2) == this%deg(cand2))then
!							conta_falhas1 = conta_falhas1 + 1
!							cycle Principal
!						endif
!
!		e nem essa
!
!						if(degAux(cand1) == this%deg(cand1)) &
!						exit Principal
!#######################################################################





						
!#######################################################################
!		Aqui evito autoconexoes e ligacoes repetidas.
!		Se nao ha ligacoes repetidas apos o 'else', ligaStubs
!		muda de valor, permitindo sair do loop .
!		Toda vez que a conexao falha, contamos para entao resetar
!		a semente, caso ocorra um numero muito grande de falhas.
!#######################################################################


						if(cand2 == cand1)then
							
							conta_falhas1 = conta_falhas1 + 1
							cycle Principal
						else
							do k = this%aux(cand1), this%aux(cand1) + degAux(cand1) - 1
								if(this%listAdj(k) == cand2) then
									conta_falhas1 = conta_falhas1 + 1
									cycle Principal
								endif
							enddo
							ligouStubs = .True.
						endif
										
						this%listAdj(this%aux(cand1)+degAux(cand1)) = cand2
						this%listAdj(this%aux(cand2)+degAux(cand2)) = cand1
						degAux(cand1) = degAux(cand1) + 1
						degAux(cand2) = degAux(cand2) + 1
						
						this%edge = this%edge + 1


!#######################################################################
!		Se escolhi gerar a lista de adjacencias pro Gephi,
!		ele gera, caso contrario, ignora.
!#######################################################################
						

						this%matriz(this%edge, 1) = cand1
						this%matriz(this%edge, 2) = cand2

!#######################################################################
!	Se stub2 vier depois de stub1, nao podemos eliminar o stub1
!	primeiro pois, se stub2 for lastC, colocaremos cand2 (stub2)
!	na posicao stub1 e, quando formos eliminar o stub2, ele estara
!	virtualmente fora da lista lC. Dito isso, ao colocarmos o
!	lastC atual na posicao stub2, estaremos tirando um provavel lastC,
!	que nao foi ligado ainda, da lista lC.
!	Situacao similar ocorre se pensarmos inicialmente no stub1.
!#######################################################################

						if(stub1 > stub2)then							
							lC(stub1) = lC(lastC)
							lastC = lastC - 1

							lC(stub2) = lC(lastC)
							lastC = lastC - 1
						else
							lC(stub2) = lC(lastC)
							lastC = lastC - 1

							lC(stub1) = lC(lastC)
							lastC = lastC - 1
						endif	
			enddo Principal
	enddo loop_stubs

!#######################################################################
!	Devido a possibilidade de o modelo UCM poder deixar stubs soltos,
!	vamos 'confirmar' a rede. Neste caso, criar a matriz de conexoes
!	passa a nao ser mais opcional, mas obrigatorio.
!
!#######################################################################			


			
!			write(*,*) "Sobraram ", nStubs - this%edge," stubs"

!#######################################################################
!	Se sobraram stubs, vamos reformular a lista de graus.
!	A lista e zerada e a partir da matriz de adjacencias,
!	nos recontamos os graus.
!#######################################################################			

			if((nStubs - this%edge) > 0)then
			
			!###########################################################
			!	Se o numero de conexoes nao eh o previsto,
			!	zeramos a lista de graus e recontamos a lista de
			!	conexoes atraves da matriz de conexoes.
			!###########################################################
				
!			write(*,*) " "
!			write(*,*) "Faremos a reestruturacao das listas da rede"
!			write(*,*) " "
!			write(*,*) "Parametros: "
!			write(*,*) " "
!			write(*,*) "Numero de nos: ", this%nodes
!			write(*,*) "Numero de ligacoes: ", this%edge
!			write(*,*) "Grau maximo: ", this%degMax
!			write(*,*) "Gama: ", gammam
!			write(*,*) "Probabilidade de ligar grau 2: ", probs
!			write(*,*) "Grau medio: ", this%degMean 
!			write(*,*) " "
		

				this%deg = 0
				
				do i = 1, this%edge
					cand1 = this%matriz(i,1)
					this%deg(cand1) = this%deg(cand1) + 1
					
					cand2 =	this%matriz(i,2)
					this%deg(cand2) = this%deg(cand2) + 1
				enddo

			!###########################################################
			!	A soma de graus e recontada
			!###########################################################			
				

!				write(*,*) " "
!				write(*,*) "Contamos ", this%sumDeg, " stubs usados em ligacoes na rede ..."
!				write(*,*) " "

			!###########################################################
			!	Guardamos o numero original de nos para testar se houve
			!	discrepancia.
			!###########################################################			
				
				N1 = this%nodes
				

			!###########################################################
			!	Zeramos o numero de nos original e recontamos a partir da lista
			!	de graus
			!###########################################################			
				this%nodes = 0
				do i = 1, N1
					if(this%deg(i) > 0) this%nodes = this%nodes + 1
				enddo

!				write(*,*) " "
!				write(*,*) "Contamos ", this%nodes, " nos ligantes ate agora..."
!				write(*,*) " "

!#######################################################################
!	Se o numero de nos decaiu, desalocamos e alocamos a lista
!	de graus e a lista auxiliar e recontamos a lista de graus.
!#######################################################################			

! Minha this%matriz ainda guarda os indices originais dos nos.
! Quando this%nodes < N1, da merda.
! Minha ideia eh reindexar os nos do de indice menor para o de indice menor.
! Mas como?


				
				if(this%nodes < N1)then

!#######################################################################			
!	Mapeamento dos indices antigos para os atuais

					if(allocated(lista_auxiliar_limpeza)) deallocate(lista_auxiliar_limpeza)
						allocate(lista_auxiliar_limpeza(N1))
						lista_auxiliar_limpeza = 0
					
					do i = 1, this%edge
						lista_auxiliar_limpeza(this%matriz(i, 1)) = 1
						lista_auxiliar_limpeza(this%matriz(i, 2)) = 1
					enddo	
					
					do i = N1, 1, -1
						if(lista_auxiliar_limpeza(i) > 0 )then
							lista_auxiliar_limpeza(i) = this%nodes
							this%nodes = this%nodes - 1
						endif  
					enddo	
					
					do i = 1, this%edge
						cand1 = this%matriz(i,1)
						this%matriz(i,1) = lista_auxiliar_limpeza(cand1)
						cand2 = this%matriz(i,2)
						this%matriz(i,2) = lista_auxiliar_limpeza(cand2)
					enddo
					
					this%nodes = maxval(this%matriz)
!#######################################################################						
							
					if(allocated(this%deg)) deallocate(this%deg)
					allocate(this%deg(this%nodes))
					this%deg = 0
						
					if(allocated(this%aux)) deallocate(this%aux)
					allocate(this%aux(this%nodes))
					
					do i = 1, this%edge
						cand1 =	this%matriz(i,1)
						this%deg(cand1) = this%deg(cand1) + 1

						cand2 =	this%matriz(i,2)
						this%deg(cand2) = this%deg(cand2) + 1
					enddo
					
				endif

!#######################################################################
!	Recontamos a lista auxiliar.
!#######################################################################			

				
				this%aux(1) = 1
					
				do i = 2, this%nodes
					this%aux(i) = this%aux(i-1) + this%deg(i-1)
				enddo

!#######################################################################
!	Finalmente, recontamos a lista de adjacencias.
!#######################################################################			

				this%sumDeg = sum(this%deg)

				if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
				this%deg = 0
				
!#######################################################################
!	Para isto, re-zeramos a lista de graus e comecamos a contar
!	tudo de novo... isso ate achar um metodo melhor.
!#######################################################################							
				do i = 1, this%edge
					
					cand1 =	this%matriz(i,1)
					cand2 =	this%matriz(i,2)
					this%listAdj(this%aux(cand1) + this%deg(cand1)) = cand2
					
					this%deg(cand1) = this%deg(cand1) + 1
					
					this%listAdj(this%aux(cand2) + this%deg(cand2)) = cand1
					
					this%deg(cand2) = this%deg(cand2) + 1					
				enddo
			endif

			this%degMax = maxval(this%deg)
			this%degMin = minval(this%deg)
			
			write(*,*) 'O grau maximo apos a ligacao da rede eh : ', this%degMax
			write(*,*) 'O grau minimo apos a ligacao da rede eh : ', this%degMin
						
			if(escreveListaConexoes)then
				open(113,file='lista_conexoes.csv', status='unknown')

				do i = 1, this%edge
					write(113,*) this%matriz(i,1), ",", this%matriz(i,2)
				enddo
				close(113)
			endif
			
			this%degMean = 1.0_dp * this%sumDeg/this%nodes
			
			write(*,*) 'Este e o grau medio apos a ligacao da rede: ', this%degMean

						
	end subroutine		


!#######################################################################
!	Subrotinas para redes RRN e misturas delas
!#######################################################################

!#######################################################################
!		Setter pra inicializar junto ki.
!		No inicializador padrao nao eh possivel.
!#######################################################################

		subroutine iniciaGrafoRRN(this, N, degRRN, seed, escreveListaConexoes)
			!-----------------------------------------------------------
			class(grafoRRN) :: this
			integer, intent(in) :: N
			integer, intent(in) :: degRRN
			integer :: seed
			logical :: escreveListaConexoes
			!-----------------------------------------------------------
			real(dp) :: harvest
			integer :: i, the_chosen
			type(rndgen) :: gerador
			!-----------------------------------------------------------
		    call iniciaGrafo(this, N)
			!-----------------------------------------------------------			
			this%degRRN = degRRN
			this%degMean = degRRN
			this%degMin = degRRN
			this%degMax = degRRN
			!-----------------------------------------------------------			
			this%aux(1) = 1
			
			do i = 2, this%nodes
				this%aux(i) = this%aux(i-1) + this%degRRN
			enddo
			
			do i = 1, this%nodes
				this%deg(i) = degRRN
			enddo
			
			this%sumDeg = sum(this%deg)


!#######################################################################
!	Se o numero total de stubs nao eh par, eu somo 1 unidade a um
!	no tomado aleatoriamente da rede.
!	Selecao aleatoria ainda nao implementada.
!#######################################################################
			
			
			if(mod(this%sumDeg, 2) > 0) then
				!##############################################
				!	Aqui uso o gerador do Fortran, porque
				!	nao consegui usar o do Wesley sem semente
				!##############################################
				write(*,*) "Algo errado aconteceu ao criar grafo RRN. Abortando codigo..."
				stop
			endif

			if(allocated(this%listAdj)) deallocate(this%listAdj)
				allocate(this%listAdj(this%sumDeg))
			
			if(allocated(this%matriz)) deallocate(this%matriz)
				allocate(this%matriz(this%sumDeg/2, 2))
!				this%matriz = 0

			!-----------------------------------------------------------
            call ligaUCM(this, seed, escreveListaConexoes)
			!-----------------------------------------------------------
		end subroutine
		
!#######################################################################
!	Mistura graus de redes regulares afim de se obter uma rede mista.
!#######################################################################

	subroutine misturaGrafos(this, p1, k1, k2)
			real(dp), intent(in) :: p1
!			real(dp), optional :: p2
			integer, intent(in) :: k1, k2
!			integer, optional :: k3
			integer :: N1, N2, N3
			integer :: i1, i2
			class(grafo), intent(inout) :: this
			

			if(p1>1) stop "Nao eh possivel misturar grafos com o valor &
			fornecido para p1"
!			if( (p1 + p2) > 1.0_dp) stop "Escolha p1 e p2 de forma que &
!			sua soma seja menor que 1"
			if(.not.allocated(this%deg)) stop "Inicie primeiro a &
			subrotina iniciaGrafo."
			
			
			N1 = int(1.0_dp * p1 * this%nodes)
!			N2 = int(1.0_dp * p2 * this%nodes)
			
			do i1 = 1, N1
				this%deg(i1) = k1
			enddo
			
			do i2 = N1+1, this%nodes
				this%deg(i2) = k2	
			enddo
			
!			N3 = 1.0_dp * (1.0_dp -(p1 + p2)) * this%nodes
			
			this%sumDeg = sum(this%deg)
			
			If(mod(this%sumDeg, 2)/=0) then
				this%deg(this%nodes) = this%deg(this%nodes) + 1
				this%sumDeg = this%sumDeg + 1
			endif
		

                        this%aux(1) = 1
                        do i1 = 2, this%nodes
                                this%aux(i1) = this%aux(i1-1) + this%deg(i1-1)
                        enddo

			if(allocated(this%listAdj)) deallocate(this%listAdj)
			allocate(this%listAdj(this%sumDeg))
			
			this%degMean = 1._dp * sum(this%deg)/this%nodes
						
		end subroutine



				
end module
