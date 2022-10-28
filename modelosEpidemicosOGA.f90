!########################################################################################################
!	Modulo que faz a dinamica
!########################################################################################################

module dynamics
	use geraRede
	use mod_rndgen
	use mod_tools	
	implicit none
	
	
	!################################################################################################
	!	Listas dinamicas
	!################################################################################################

	integer, allocatable :: infList(:), infList1(:), infList2(:)		!, susList(:), infSusList(:)
!	integer(kind=1), allocatable :: sigma(:)
	
	!################################################################################################
	!	Variaveis dinamicas globais
	!################################################################################################

	
	real(dp), private, parameter :: mu = 1._dp
	
	real(dp) :: lambda
	real(dp) :: dlambda, lambda_i, lambda_f
	

	!################################################################################################
	!	Para inicializar
	!################################################################################################
	
	real(dp) :: nInf0
	
	!################################################################################################
	!	Variaveis dinamicas
	!################################################################################################
	
	integer :: 	nInf1, nInf2, nInf, nInf_n_redundante,&
	 nInfMax, nInfMax_n_redundante, nInf_redundante, nInf_redundante1, nInf_redundante2,&
	 nInfMax_redundante, nInfMax_redundante1, nInfMax_redundante2, nSus, nInfSus

	
	!################################################################################################
	!	Taxas de eventos
	!################################################################################################

	real(dp) :: rateTotal, lambdaTotal, muTotal

	
	!################################################################################################
	!	Probabilidades associadas aos eventos independentes
	!################################################################################################
	
	real(dp) :: m, l
	

	!################################################################################################
	!	Variaveis temporais
	!################################################################################################

	real(dp) :: t, dt, life_spam_time, life_spam_time_n_redundante
	
	integer	:: t_pos, t_Media, t_relax, tMax, t_MaxVis
	real(dp) :: t_LS, t_LS_n_redundante, t_LS_redundante, t_LS_redundante1, t_LS_redundante2
		
	real(dp), allocatable :: t_Medio(:)


	!################################################################################################
	!	Variaveis temporais auxiliares
	!################################################################################################

	
	integer, allocatable :: t_samp(:), samp_surv(:)
	integer (kind=1), allocatable :: vezes_update(:)

	!################################################################################################
	!	Variaveis epidemicas medias
	!################################################################################################

	
	real(dp), allocatable :: rho_medio(:)
	real(dp) :: rho_medioQS, rho_medioQS_n_redundante,rho_medioQS_redundante,rho_medioQS_redundante1,rho_medioQS_redundante2, &
	 rho2_medioQS, rho2_medioQS_n_redundante,rho2_medioQS_redundante, rho2_medioQS_redundante1, rho2_medioQS_redundante2, &
	 dev_rho_medioQS, dev_rho_medioQS_n_redundante, dev_rho_medioQS_redundante, dev_rho_medioQS_redundante1, dev_rho_medioQS_redundante2 
	real(dp) :: Xi, Xi_n_redundante, Xi_redundante, Xi_redundante1, Xi_redundante2 

	real(dp) :: S_Schanon, S_Schanon_n_redundante, S_Schanon_redundante, S_Schanon_redundante1, S_Schanon_redundante2
	!################################################################################################
	!	gerador de numeros aleatorios
	!################################################################################################
	
	type(rndgen) :: gen

	
	!################################################################################################
	!	variavel aleatoria
	!################################################################################################
	
	real(dp) :: prob
	
	
	!################################################################################################
	!	Variaveis auxiliares
	!################################################################################################
	integer, allocatable :: listAux(:)
	integer, allocatable :: maskara(:)
	integer :: nSamples, npts, t_up, n_up, n_upi, comp_gigante_dummy, nInf_dummy
	integer ::  n_K, nInf_pfoto
	logical :: recuperaDinamica
	logical :: imuniza_tubos
	real(dp) :: esquenta
	integer, allocatable :: infectaveis_cg(:)
    integer :: nInf1_dummy, nInf2_dummy
	integer :: W, W1, W2, W_dummy, W1_dummy, W2_dummy
        !################################################################################################
        !       Auxiliar do ioga
        
	integer, allocatable :: abaixo_media(:), acima_media(:)
	integer :: grauMax1, grauMax2
        
	!################################################################################################
	!	Vetor que vai guardar a Probabilidade quasi-estacionaria.
	!################################################################################################

	real(dp), allocatable :: Pn_QS(:), Pn_QS_n_redundante(:), Pn_QS_redundante(:),Pn_QS_redundante1(:), Pn_QS_redundante2(:) 
	real(dp) :: sumPQS
	
	!##########################################	
	!	Mudanda de obj pra procedural
		integer:: tamRede, grauMax, grauMin
		integer :: grauMin2
		integer, allocatable :: lista_aux(:), lista_adj(:), grau(:)
	!##########################################	
	
	!###################################################################
	!	Listas para o Metodo Quase-Estacionario Padrao
	!###################################################################
	
		integer, allocatable :: fotos_din(:,:), nInf_fotos_din(:)
	
	
	
	contains

!####################################################################################	
!		Modelo S-I-S e processo de contato
!####################################################################################


		subroutine allocaMedias_An_Temporal(this, tmax)
			class(grafo), intent(in) :: this
			integer, intent(in) :: tmax
			
			if(allocated(rho_medio)) deallocate(rho_medio)
			allocate(rho_medio(tMax))
			rho_medio = 0._dp

			
			!################################################################################
			!	Lista para guardar medias dos tempos de eventos
			!################################################################################
			
			
			if(allocated(t_Medio)) deallocate(t_Medio)
			allocate(t_Medio(tMax))
			t_medio = 0._dp
			
			if(allocated(t_samp)) deallocate(t_samp)										! Vai guardar o numero de amostras
			allocate(t_samp(tMax))															! que participam da dinamica num dado tempo discreto.
			t_samp = 0
			
			if(allocated(samp_surv)) deallocate(samp_surv)									
			allocate(samp_surv(tMax))
			samp_surv = 0
			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(this%nodes))
				

			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))
									
			t_MaxVis = 0
		end subroutine


!####################################################################################

		subroutine allocaMedias_An_imun_tubos_e_folhas(this)
			class(grafo) :: this
!			integer :: comp_gigante_dummy, nInf_dummy
			integer :: i19, i20, i21, n_infectaveis
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################


			n_infectaveis=0
			
l1:			do i20 = 1, this%nodes
	
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l1
				if(this%deg(i20) > 2) n_infectaveis = n_infectaveis + 1
	
			enddo l1
			
			if(allocated(infectaveis_cg)) deallocate(infectaveis_cg)
				allocate(infectaveis_cg(n_infectaveis))
				
			n_infectaveis = 0

!#######################################################################
!	O número de infectados sera calculado com exatidao
!#######################################################################
						
			nInfMax_n_redundante = 0
			
l2:			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l2
					if(this%deg(i20) > 2)then
						n_infectaveis = n_infectaveis + 1
						infectaveis_cg(n_infectaveis) = i20
!						if(this%deg(i20) == 2) nInfMax_redundante2 = nInfMax_redundante2 + 1
						if(this%deg(i20) > 2) nInfMax_n_redundante = nInfMax_n_redundante + 1							
					endif
			enddo l2


!#######################################################################			
			nInfMax = n_infectaveis
!#######################################################################
			
			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
	
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
			
			if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
			allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						

!			if(allocated(Pn_QS_redundante2)) deallocate(Pn_QS_redundante2)
!			allocate(Pn_QS_redundante2(0:nInfMax_redundante2))									

				!#######################################################
				!	Pode sobrar alguem de grau 2 na decomposicao
				!	k-onion
				!#######################################################			
	
											
		end subroutine


!####################################################################################

		subroutine condicaoInicial_CG_imun_tubos_e_folhas(this)
			class(grafo) :: this
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: i1, j1, j2, k1, lastCand
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			
			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

						sigma = 0																			! Todos os nos sao suscetiveis				
						
						infList = 0

		!#######################################################################################									

						k1 = 1
eh_esse:				do i1 = 1, this%nodes

!#######################################################################
!	Imuniza quem não é da componente gigante
!#######################################################################
							
							if(lista_de_clusters(i1) /= i_comp_gigante)then
								sigma(i1) = 2
								cycle eh_esse
							endif	

							if(this%deg(i1) < 3)then
								sigma(i1) = 2
								cycle eh_esse
							endif
										
							
							infList(k1) = i1
							sigma(i1) = 1
							k1 = k1 + 1
						enddo eh_esse
						
						nInf = nInfMax
						
						nInf_n_redundante = nInfMax_n_redundante
						
!						nInf_redundante2 = nInfMax_redundante2
						
						t = 0._dp;
						
						t_pos = 1;
						
		end subroutine	

!##


!####################################################################################
!	Imuniza nos de grau dois dos tubos, mas que sejam vizinhos apenas de grau 2.
!####################################################################################

		subroutine allocaMedias_An_imun_so_folhas(this)
			class(grafo) :: this
!			integer :: comp_gigante_dummy, nInf_dummy
			integer :: i19, i20, i21, n_infectaveis
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################


			n_infectaveis=0
			
l1:			do i20 = 1, this%nodes
	
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l1
				if(this%deg(i20) > 1) n_infectaveis = n_infectaveis + 1
	
			enddo l1
			
			if(allocated(infectaveis_cg)) deallocate(infectaveis_cg)
				allocate(infectaveis_cg(n_infectaveis))
				
			n_infectaveis = 0

!#######################################################################
!	O número de infectados sera calculado com exatidao
!#######################################################################
			nInfMax_redundante2 = 0
			nInfMax_n_redundante = 0
			
l2:			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l2
					if(this%deg(i20) > 1)then
						n_infectaveis = n_infectaveis + 1
						infectaveis_cg(n_infectaveis) = i20
						if(this%deg(i20) == 2) nInfMax_redundante2 = nInfMax_redundante2 + 1
						if(this%deg(i20) > 2) nInfMax_n_redundante = nInfMax_n_redundante + 1							
					endif
			enddo l2


!#######################################################################			
			nInfMax = n_infectaveis
!#######################################################################
			
			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
	
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
			
			if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
			allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						

			if(allocated(Pn_QS_redundante2)) deallocate(Pn_QS_redundante2)
			allocate(Pn_QS_redundante2(0:nInfMax_redundante2))									

				!#######################################################
				!	Pode sobrar alguem de grau 2 na decomposicao
				!	k-onion
				!#######################################################			
	
											
		end subroutine

!####################################################################################
		subroutine condicaoInicial_CG_imun_so_folhas(this)
			class(grafo) :: this
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: i1, j1, j2, k1, lastCand
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			
			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

						sigma = 0																			! Todos os nos sao suscetiveis				
						
						infList = 0

		!#######################################################################################									

						k1 = 1
eh_esse:				do i1 = 1, this%nodes

!#######################################################################
!	Imuniza quem não é da componente gigante
!#######################################################################
							
							if(lista_de_clusters(i1) /= i_comp_gigante)then
								sigma(i1) = 2
								cycle eh_esse
							endif	

							if(this%deg(i1) == 1)then
								sigma(i1) = 2
								cycle eh_esse
							endif
							
							infList(k1) = i1
							sigma(i1) = 1
							k1 = k1 + 1
						enddo eh_esse
						
						nInf = nInfMax
						
						nInf_n_redundante = nInfMax_n_redundante
						
						nInf_redundante2 = nInfMax_redundante2
						
						t = 0._dp;
						
						t_pos = 1;
						
		end subroutine	

!####################################################################################

		subroutine allocaMedias_An_imun_so_tubos(this)
			class(grafo) :: this
!			integer :: comp_gigante_dummy, nInf_dummy
			integer :: i19, i20, i21, n_infectaveis
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################


			n_infectaveis=0
l1:			do i20 = 1, this%nodes
				
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l1
				if(this%deg(i20) /= 2) n_infectaveis = n_infectaveis + 1
					
			enddo l1
			
			if(allocated(infectaveis_cg)) deallocate(infectaveis_cg)
				allocate(infectaveis_cg(n_infectaveis))
				


!#######################################################################
!	O número de infectados sera calculado com exatidao
!#######################################################################
			n_infectaveis = 0
			nInfMax_redundante1 = 0
			nInfMax_n_redundante = 0
!			nInfMax_redundante2 = 0			

l2:			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle l2
					if(this%deg(i20) /=2)then
						
						n_infectaveis = n_infectaveis + 1
						
						infectaveis_cg(n_infectaveis) = i20
						
						if(this%deg(i20) == 1)then
							nInfMax_redundante1 = nInfMax_redundante1 + 1
						elseif(this%deg(i20) > 2)then
							nInfMax_n_redundante = nInfMax_n_redundante + 1
						endif
														
					endif
			enddo l2


!#######################################################################			
			nInfMax = n_infectaveis
!#######################################################################
			
			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
							
			if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
			allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						

			if(allocated(Pn_QS_redundante1)) deallocate(Pn_QS_redundante1)
			allocate(Pn_QS_redundante1(0:nInfMax_redundante1))									

				!#######################################################
				!	Pode sobrar alguem de grau 2 na decomposicao
				!	k-onion
				!#######################################################
																			
		end subroutine
			

!#######################################################################

		subroutine condicaoInicial_CG_imun_so_tubos(this)
			class(grafo) :: this
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: i1, j1, j2, k1, lastCand
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			
			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

						sigma = 0																			! Todos os nos sao suscetiveis				
			
						infList = 0
						
		!#######################################################################################									

						k1 = 1
eh_esse:				do i1 = 1, this%nodes

!#######################################################################
!	Imuniza quem não é da componente gigante
!#######################################################################
							
							if(lista_de_clusters(i1) /= i_comp_gigante)then
								sigma(i1) = 2
								cycle eh_esse
							endif	

							if(this%deg(i1) == 2)then
								sigma(i1) = 2
								cycle eh_esse
							endif
							
							infList(k1) = i1
							sigma(i1) = 1
							k1 = k1 + 1
						enddo eh_esse
						nInf = nInfMax
						nInf_redundante1 = nInfMax_redundante1
						nInf_n_redundante = nInfMax_n_redundante
						t = 0._dp;
						t_pos = 1;
		end subroutine	


!#######################################################################
				!temporariamente aqui
			
			
!#######################################################################		
		subroutine condicaoInicial_CG1(this)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			class(grafo) :: this
			integer :: i1, j1, k1, lastCand, nInf_dummy
!			logical, intent(in) :: imuniza_tubos
			

		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

				sigma = 0																			! Todos os nos sao suscetiveis
																								
				k1 = 0

lupy:			do i1 = 1, this%nodes

				if(lista_de_clusters(i1) /= i_comp_gigante)then
					sigma(i1) = 2
					cycle lupy
				endif	
				
				k1 = k1 + 1
				infList(k1) = i1
				sigma(i1) = 1
				
			enddo lupy
						
			nInf = nInfMax

			nInf_n_redundante = nInfMax_n_redundante
			nInf_redundante1 = nInfMax_redundante1
			nInf_redundante2 = nInfMax_redundante2

			t = 0.0_dp;
			t_pos = 1;

		end subroutine
				

!####################################################################################
!		Processo SIS com condicao de contorno reflexiva
!####################################################################################


		subroutine sisProcessRBS_imun_folhas(this, com_comp_gigante, seed)


!###################################################################################
!	A rede, ou substrato
!###################################################################################

!#######################################################################
			class(grafo), intent(in) :: this
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
			logical, intent(in) :: com_comp_gigante
!#######################################################################
				 
				 n_K = 0

				 do i1 = 1, nInf
					n_K = n_K + this%deg(infList(i1))
!					n_K = n_K + grau2(infList(i1))
				 enddo
				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp
				Pn_QS_redundante2 = 0.0_dp				
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		


			
loopDinamico:	do while(t <= tMax)
				
				
				
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = 1.0_dp * nInf * mu
				
				lambdaTotal = 1.0_dp * n_K * lambda

				rateTotal = muTotal + lambdaTotal

				m = 1.0_dp * muTotal/rateTotal
						
	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt

	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################	

				
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante2(nInf_redundante2) = Pn_QS_redundante2(nInf_redundante2) + dt										
				endif

				
	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	


					if(nInf == 1)then
						sigma(infList(1)) = 0
						n_K = 0					
						
!						if(grau2(infList(1)) > 2) nInf_n_redundante = 0
!						if(grau2(infList(1)) == 2) nInf_redundante2 = 0						

						if(this%deg(infList(1)) > 2) nInf_n_redundante = 0
						if(this%deg(infList(1)) == 2) nInf_redundante2 = 0						
																				
l_heal_last:            do
							i1 = gen%int(1,size(infectaveis_cg))
							if(sigma(infectaveis_cg(i1)) == 0) exit l_heal_last
						enddo l_heal_last
						
						infList(1) = infectaveis_cg(i1)			
						
!						if(grau2(infList(1)) > 2) nInf_n_redundante = 1
!						if(grau2(infList(1)) == 2) nInf_redundante2 = 1						

						if(this%deg(infList(1)) > 2) nInf_n_redundante = 1
						if(this%deg(infList(1)) == 2) nInf_redundante2 = 1						
						
						sigma(infList(1)) = 1
						
!						n_K = this%deg(infList(1))							!n_K + grau(infList(1))
						n_K = this%deg(infList(1))							!n_K + grau(infList(1))												
						cycle loopDinamico
					endif

	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
					
					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
!					n_K = n_K - grau2(infList(i1))
										
					if(this%deg(infList(i1)) > 2) nInf_n_redundante = nInf_n_redundante - 1					
					if(this%deg(infList(i1)) == 2) nInf_redundante2 = nInf_redundante2 - 1

					infList(i1) = infList(nInf)
					nInf = nInf - 1

	!###################################################################################################################
	!	Como soh ha dois eventos possiveis, cura e infeccao, no caso de a cura nao ser sorteada, sobra a infeccao
	!	automaticamente.
	!	Um no infectado eh sorteado, estabelecemos a probabilidade dele infectar baseado no seu grau de infeccao
	!	probInf = 1.0_dp * lambda * this%deg(infList(i1))/ 
	!###################################################################################################################	


					
		   else	
loop_infecta:		do																	! Se o processo da vez eh infeccao.
					i1 = gen%int(1, nInf)
					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle loop_infecta

					j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

					if(sigma(this%listAdj(j1)) == 0)then
						sigma(this%listAdj(j1)) = 1

						n_K = n_K + this%deg(this%listAdj(j1))
!						n_K = n_K + grau2(this%listAdj(j1))
						nInf = nInf + 1			
						infList(nInf) = this%listAdj(j1)
						
						if(this%deg(infList(nInf)) > 2) nInf_n_redundante = nInf_n_redundante + 1
						if(this%deg(infList(nInf)) == 2) nInf_redundante2 = nInf_redundante2 + 1												
												
					endif						
					exit loop_infecta
			enddo loop_infecta
		   endif		
		   					
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para Pn_QS quando para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################				
								
			Pn_QS = 1.0_dp * Pn_QS / sumPQS
						
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp
			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo
			
			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante
			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((1.0_dp * comp_gigante)**2.0_dp)
			
			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS

!#######################################################################
! Para a parte nao redundante

			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
						
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS_n_redundante(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp
			
			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			
!			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante_din
			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante
!			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante_din)**2.0_dp)
			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((1.0_dp * nInfMax_n_redundante)**2.0_dp)
						
			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
!			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante				
			
			!###########################################################
			!	Nos de grau 2
			
			
			Pn_QS_redundante2 = 1.0_dp * Pn_QS_redundante2 / sumPQS
						
			if(Pn_QS_redundante2(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS_redundante2(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS_redundante2 = 0.0_dp
			rho2_medioQS_redundante2 = 0.0_dp
			S_Schanon_redundante2 = 0.0_dp
			
			do i1 = 1, nInfMax_redundante2
				if(Pn_QS_redundante2(i1) > 0.0_dp)then
					rho_medioQS_redundante2 = 1.0_dp * i1 * Pn_QS_redundante2(i1) + rho_medioQS_redundante2						! n medio
					rho2_medioQS_redundante2 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante2(i1) + rho2_medioQS_redundante2
					S_Schanon_redundante2 = S_Schanon_redundante2 - Pn_QS_redundante2(i1) * log(Pn_QS_redundante2(i1))
				endif
			enddo
			
!			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante_din
			rho_medioQS_redundante2 = 1.0_dp * rho_medioQS_redundante2 / nInfMax_redundante2
!			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante_din)**2.0_dp)
			rho2_medioQS_redundante2 = 1.0_dp * rho2_medioQS_redundante2/((1.0_dp * nInfMax_redundante2)**2.0_dp)
						
			dev_rho_medioQS_redundante2 = (rho2_medioQS_redundante2 - rho_medioQS_redundante2**2.0_dp)**0.5_dp
			
!			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
			Xi_redundante2 = 1.0_dp * nInfMax_redundante2 * (rho2_medioQS_redundante2 - (rho_medioQS_redundante2**2.0_dp))/rho_medioQS_redundante2				
				
		end subroutine		


!####################################################################################
		subroutine sisProcessRBS_imun_tubos(this, com_comp_gigante, seed)


!###################################################################################
!	A rede, ou substrato
!###################################################################################

!#######################################################################
			class(grafo), intent(in) :: this
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
			logical, intent(in) :: com_comp_gigante
!#######################################################################
				 
				 n_K = 0

				 do i1 = 1, nInf
					n_K = n_K + this%deg(infList(i1))
				 enddo

!				 do i1 = 1, nInf
!					n_K = n_K + grau2(infList(i1))
!				 enddo

				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp
				Pn_QS_redundante1 = 0.0_dp
				
							
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		


			
loopDinamico:	do while(t <= tMax)
				
				
				
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = (1.0_dp * nInf) * mu

				lambdaTotal = (1.0_dp * n_K) * lambda

				rateTotal = muTotal + lambdaTotal

				m = 1.0_dp * muTotal/rateTotal

				
				
	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt

	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################	

										
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante1(nInf_redundante1) = Pn_QS_redundante1(nInf_redundante1) + dt					
				endif

				
	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	


					if(nInf == 1)then
						sigma(infList(1)) = 0
						n_K = 0					

!						if(grau2(infList(1)) > 2) nInf_n_redundante = 0
!						if(grau2(infList(1)) == 1) nInf_redundante1 = 0						

						if(this%deg(infList(1)) > 2)then
							nInf_n_redundante = 0
						elseif(this%deg(infList(1)) == 1)then
							nInf_redundante1 = 0
						endif						

l_heal_last:            do
							i1 = gen%int(1,size(infectaveis_cg))
							if(sigma(infectaveis_cg(i1)) == 0) exit l_heal_last
						enddo l_heal_last
						
						infList(1) = infectaveis_cg(i1)			

						if(this%deg(infList(1)) > 2)then
							nInf_n_redundante = 1
						elseif(this%deg(infList(1)) == 1)then
							nInf_redundante1 = 1						
						endif
						
!						if(grau2(infList(1)) > 2) nInf_n_redundante = 1
!						if(grau2(infList(1)) == 1) nInf_redundante1 = 1						
								
						sigma(infList(1)) = 1

						n_K = this%deg(infList(1))
!						n_K = grau2(infList(1))												

						cycle loopDinamico
					endif

	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
					
					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0

					n_K = n_K - this%deg(infList(i1))

!					n_K = n_K - grau2(infList(i1))					

					if(this%deg(infList(i1)) > 2) then
						nInf_n_redundante = nInf_n_redundante - 1
					elseif(this%deg(infList(i1)) == 1) then
						nInf_redundante1 = nInf_redundante1 - 1					
					endif

					infList(i1) = infList(nInf)
																					
					nInf = nInf - 1
					
	!###################################################################################################################
	!	Como soh ha dois eventos possiveis, cura e infeccao, no caso de a cura nao ser sorteada, sobra a infeccao
	!	automaticamente.
	!	Um no infectado eh sorteado, estabelecemos a probabilidade dele infectar baseado no seu grau de infeccao
	!	probInf = 1.0_dp * lambda * this%deg(infList(i1))/ 
	!###################################################################################################################	


					
		   else	
loop_infecta:		do																	! Se o processo da vez eh infeccao.

					i1 = gen%int(1, nInf)
					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle loop_infecta

					j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

					if(sigma(this%listAdj(j1)) == 0)then

						sigma(this%listAdj(j1)) = 1

						nInf = nInf + 1					
						infList(nInf) = this%listAdj(j1)	

						n_K = n_K + this%deg(infList(nInf))
!						n_K = n_K + grau2(this%listAdj(j1))
						
												
							
!						if(grau2(infList(1)) > 2) nInf_n_redundante = nInf_n_redundante + 1
!						if(grau2(infList(1)) == 1) nInf_redundante1 = nInf_redundante1	+ 1					

						if(this%deg(infList(nInf)) > 2) then
							nInf_n_redundante = nInf_n_redundante + 1
						elseif(this%deg(infList(nInf)) == 1)then
							nInf_redundante1 = nInf_redundante1 + 1					
						endif
						
					endif						
					
					exit loop_infecta
			enddo loop_infecta
		   endif		
		   					
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para Pn_QS quando para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################				
			! Para a rede toda
			Pn_QS = 1.0_dp * Pn_QS / sumPQS
						
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp
			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo
			
!			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante_din
			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante
!			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante_din)**2.0_dp)
			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((1.0_dp * comp_gigante)**2.0_dp)
						
			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
!			Xi = 1.0_dp * comp_gigante_din * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS				

!#######################################################################
! Para a parte nao redundante

			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
						
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS = 1.0_dp/Pn_QS_n_redundante(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp
			
			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			
!			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante_din
			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante
!			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante_din)**2.0_dp)
			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((1.0_dp * nInfMax_n_redundante)**2.0_dp)
						
			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
!			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante				
			
			!###########################################################
			!	Nos de grau 1
			
			
			Pn_QS_redundante1 = 1.0_dp * Pn_QS_redundante1 / sumPQS
						
			if(Pn_QS_redundante1(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS_redundante1(1)
			else
				t_LS = tMax
			endif
			
			rho_medioQS_redundante1 = 0.0_dp
			rho2_medioQS_redundante1 = 0.0_dp
			S_Schanon_redundante1 = 0.0_dp
			
			do i1 = 1, nInfMax_redundante1
				if(Pn_QS_redundante1(i1) > 0.0_dp)then
					rho_medioQS_redundante1 = 1.0_dp * i1 * Pn_QS_redundante1(i1) + rho_medioQS_redundante1						! n medio
					rho2_medioQS_redundante1 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante1(i1) + rho2_medioQS_redundante1
					S_Schanon_redundante1 = S_Schanon_redundante1 - Pn_QS_redundante1(i1) * log(Pn_QS_redundante1(i1))
				endif
			enddo
			
!			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante_din
			rho_medioQS_redundante1 = 1.0_dp * rho_medioQS_redundante1 / nInfMax_redundante1
!			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante_din)**2.0_dp)
			rho2_medioQS_redundante1 = 1.0_dp * rho2_medioQS_redundante1/ ((1.0_dp * nInfMax_redundante1)**2.0_dp)
						
			dev_rho_medioQS_redundante1 = (rho2_medioQS_redundante1 - rho_medioQS_redundante1**2.0_dp)**0.5_dp
			
!			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
			Xi_redundante1 = 1.0_dp * nInfMax_redundante1 * (rho2_medioQS_redundante1 - (rho_medioQS_redundante1**2.0_dp))/rho_medioQS_redundante1				
			
			
		end subroutine		
				
				
!#######################################################################
				
		subroutine allocaMedias_An_QS(this)
			class(grafo) :: this

			integer :: i20
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################
			
			
			!###########################################################
			!	
				nInfMax = comp_gigante
			!###########################################################			

			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
			
			nInfMax_redundante = 0
			nInfMax_n_redundante = 0

!#######################################################################

			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle				
				if(this%deg(i20) < 3) nInfMax_redundante = nInfMax_redundante + 1
				if(this%deg(i20) >= 3) nInfMax_n_redundante = nInfMax_n_redundante + 1

			enddo
						
			if(nInfMax_n_redundante > 0 )then
				if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
				allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						
			endif

			if(nInfMax_redundante > 0)then
				if(allocated(Pn_QS_redundante)) deallocate(Pn_QS_redundante)
				allocate(Pn_QS_redundante(0:nInfMax_redundante))									
			endif
			
		end subroutine

!#######################################################################
				
		subroutine allocaMedias_An_QS2(this)
			class(grafo) :: this

			integer :: i20
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################
			
			
			!###########################################################
			!	
				nInfMax = comp_gigante
			!###########################################################			

			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
			
			nInfMax_redundante1 = 0
			nInfMax_redundante2 = 0			
			nInfMax_n_redundante = 0

!#######################################################################

			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle				
				if(this%deg(i20) == 1) nInfMax_redundante1 = nInfMax_redundante1 + 1
				if(this%deg(i20) == 2) nInfMax_redundante2 = nInfMax_redundante2 + 1				
				if(this%deg(i20) > 2) nInfMax_n_redundante = nInfMax_n_redundante + 1
			enddo
						
				if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
				allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						

				if(allocated(Pn_QS_redundante1)) deallocate(Pn_QS_redundante1)
				allocate(Pn_QS_redundante1(0:nInfMax_redundante1))									

				if(allocated(Pn_QS_redundante2)) deallocate(Pn_QS_redundante2)
				allocate(Pn_QS_redundante2(0:nInfMax_redundante2))									


			
		end subroutine

!#######################################################################
!	Aloca tambem listas importantes ao metodo QS Padrao

		subroutine allocaMedias_An_QS_Padrao(this, nfotos)
			class(grafo) :: this
			integer, intent(in) :: nfotos
			integer :: i20
			
			
			!################################################################################
			!	Listas de infectados, estados e distribuicao de probabilidades
			!################################################################################
			
			
			!###########################################################
			!	
				nInfMax = comp_gigante
			!###########################################################			

			if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."			
			if(allocated(infList)) deallocate(infList)
			allocate(infList(nInfMax))
						
			if(allocated(sigma)) deallocate(sigma)
			allocate(sigma(this%nodes))		
			
			if(allocated(Pn_QS)) deallocate(Pn_QS)
			allocate(Pn_QS(nInfMax))
			
			nInfMax_redundante1 = 0
			nInfMax_redundante2 = 0			
			nInfMax_n_redundante = 0
			
!			nInf_pfoto = this%nodes/5
			nInf_pfoto = this%nodes

			!###########################################################
			if(allocated(fotos_din)) deallocate(fotos_din)
					allocate(fotos_din(nInf_pfoto, nfotos))

			if(allocated(nInf_fotos_din)) deallocate(nInf_fotos_din)
					allocate(nInf_fotos_din(nfotos))					
			!###########################################################
			
!#######################################################################

			do i20 = 1, this%nodes
				if(lista_de_clusters(i20) /= i_comp_gigante) cycle				
				if(this%deg(i20) == 1) nInfMax_redundante1 = nInfMax_redundante1 + 1
				if(this%deg(i20) == 2) nInfMax_redundante2 = nInfMax_redundante2 + 1				
				if(this%deg(i20) >  2) nInfMax_n_redundante = nInfMax_n_redundante + 1
			enddo
						
			if(allocated(Pn_QS_n_redundante)) deallocate(Pn_QS_n_redundante)
			allocate(Pn_QS_n_redundante(0:nInfMax_n_redundante))						

			if(allocated(Pn_QS_redundante1)) deallocate(Pn_QS_redundante1)
			allocate(Pn_QS_redundante1(0:nInfMax_redundante1))									

			if(allocated(Pn_QS_redundante2)) deallocate(Pn_QS_redundante2)
			allocate(Pn_QS_redundante2(0:nInfMax_redundante2))									
			
		end subroutine

!#######################################################################
		subroutine sisProcess_QS_Padrao(this,com_comp_gigante,nfotos, seed)

!#######################################################################
			class(grafo), intent(in) :: this
			integer, intent(in) :: nfotos
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
            real(dp), parameter :: p_deco = 0.02_dp
            real(dp) :: prob1, p_act, p_fill
            integer :: ultima_foto, n_vezes
			logical, intent(in) :: com_comp_gigante
!#######################################################################

			if( .not. allocated(fotos_din)) stop "Provavelmente a rotina allocaMedias_An_QS_Padrao nao foi executada. Faca isso."

!#######################################################################
!	nInfMax foi settado na rotina condicaoIncial_CG
!#######################################################################				 
			 if(lambda == 0.0_dp) lambda = dlambda
						
			 n_K = 0
			 
			 do i1 = 1, nInf
				n_K = n_K + this%deg(infList(i1))
			 enddo
				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp 
				Pn_QS_redundante1 = 0.0_dp 
				Pn_QS_redundante2 = 0.0_dp 
				
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		
				
				ultima_foto = 0
			
loopDinamico:	do while(t <= tMax)
							
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = 1.0_dp * nInf * mu
				lambdaTotal = 1.0_dp * n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				

	!###################################################################################################################!
	!	Precisamos deixar o sistema relaxar um pouco antes de copiar fotos. Vou usar uma janela de 1 unidade de tempo.	!
	!	Sugestao do Guilherme.																							!
	!###################################################################################################################!	
				
!				if(nInf == 0) write(*,*) "Estado absorvente alcancado. Tem nem como curar. Instante de tempo ", t


	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt
	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante1(nInf_redundante1) = Pn_QS_redundante1(nInf_redundante1) + dt
					Pn_QS_redundante2(nInf_redundante2) = Pn_QS_redundante2(nInf_redundante2) + dt					
				endif
				
				if(t > 2_dp)then
					if(ultima_foto < nfotos)then												
							prob1 = gen%rnd()
							p_fill = 1_dp * dt
							if(prob1 <= p_fill)then	
								ultima_foto = ultima_foto + 1
								nInf_fotos_din(ultima_foto) = nInf
								do i1 = 1, nInf
									fotos_din(i1, ultima_foto) = infList(i1)
								enddo
							endif	
					else							
							prob1 = gen%rnd()
							p_act = 1.0_dp * p_deco * dt						
							if(prob1 <= p_act)then
								j1 = gen%int(1, ultima_foto)
								nInf_fotos_din(j1) = nInf
								do i1 = 1, nInf
									fotos_din(i1, j1) = infList(i1)
								enddo
							endif
					endif									
				endif

	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	
			
	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
	
	
	!###################################################################
	!	Recomendacao do professor Silvio
	!###################################################################					
					if(t < 20.0_dp)then
						if(nInf == 1)then
!							write(*,*) "Ta na lanterna, no instante de tempo ", t
							sigma(infList(nInf)) = 0

							n_K = 0						
						
							nInf_n_redundante = 0; nInf_redundante1 = 0; nInf_redundante2 = 0																					
						
							n_vezes = 0
							
Lanterna:					do
								n_vezes = n_vezes + 1
								i1 = gen%int(1, this%nodes)
								if(lista_de_clusters(i1) /= i_comp_gigante) cycle Lanterna
								exit Lanterna
							enddo Lanterna

!							write(*,*) "Foi rejeitado ", n_vezes, " vezes"
							
							infList(nInf) = i1
							sigma(i1) = 1
							n_K = this%deg(i1)
							
							if(this%deg(i1) >  2) nInf_n_redundante = 1
							if(this%deg(i1) == 1) nInf_redundante1 = 1
							if(this%deg(i1) == 2) nInf_redundante2 = 1																						
							
							cycle loopDinamico			
						endif
					endif
					
!#######################################################################				
!	Se da pra curar, nos cura
!#######################################################################

					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					
					if(this%deg(infList(i1)) > 2) nInf_n_redundante = nInf_n_redundante - 1
					if(this%deg(infList(i1)) == 1) nInf_redundante1 = nInf_redundante1 - 1
					if(this%deg(infList(i1)) == 2) nInf_redundante2 = nInf_redundante2 - 1															
					
					infList(i1) = infList(nInf)
					nInf = nInf - 1

!#######################################################################################################################
!		Se o estado absorvente eh alcancado, escolhemos uma configuracao ja visitada

					if(nInf == 0)then
												
!						write(*,*) "Foi rebaixado no instante de tempo ", t
						j1 = gen%int(1,ultima_foto)
						
						nInf = nInf_fotos_din(j1)
											
						do i1 = 1, nInf
							infList(i1) = fotos_din(i1, j1)
							sigma(infList(i1)) = 1	

							if(this%deg(infList(i1)) > 2) nInf_n_redundante = nInf_n_redundante + 1
							if(this%deg(infList(i1)) == 1) nInf_redundante1 = nInf_redundante1 + 1
							if(this%deg(infList(i1)) == 2) nInf_redundante2 = nInf_redundante2 + 1														

							n_K = n_K + this%deg(infList(i1))
						enddo						
!						cycle loopDinamico
					endif
!#######################################################################################################################		
		   else	
loop_infecta:		do	

					i1 = gen%int(1, nInf)

					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle

						j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

						if(sigma(this%listAdj(j1)) == 0)then
						
							sigma(this%listAdj(j1)) = 1
							n_K = n_K + this%deg(this%listAdj(j1))
							
							nInf = nInf + 1																			! Numero de infectados aumenta.
							infList(nInf) = this%listAdj(j1)
											
							if(this%deg(this%listAdj(j1)) > 2) nInf_n_redundante = nInf_n_redundante + 1													
							if(this%deg(this%listAdj(j1)) == 1) nInf_redundante1 = nInf_redundante1 + 1
							if(this%deg(this%listAdj(j1)) == 2) nInf_redundante2 = nInf_redundante2 + 1
																											
						endif
						exit loop_infecta
			enddo loop_infecta
		   endif		
		   			
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para a Pn_QS, quanto para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################					
							
			do i1 = 1, size(Pn_QS)			
				Pn_QS(i1) = 1.0_dp * Pn_QS(i1) / sumPQS
			enddo
									
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp			
			
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif

			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo

			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante

			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((1.0_dp * comp_gigante)**2.0_dp)

			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
						


!#######################################################################
!	Estatistica nos nao redundantes
!#######################################################################
			
			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp		
			
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS_n_redundante = 1.0_dp /Pn_QS_n_redundante(1)
			else
				t_LS_n_redundante = tMax
			endif

			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			

			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante

			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((1.0_dp * nInfMax_n_redundante)**2.0_dp)

			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante

!#######################################################################
!	Estatistica nos redundantes1
!#######################################################################
			
			Pn_QS_redundante1 = 1.0_dp * Pn_QS_redundante1 / sumPQS
			rho_medioQS_redundante1 = 0.0_dp
			rho2_medioQS_redundante1 = 0.0_dp
			S_Schanon_redundante1 = 0.0_dp		

			if(nInfMax_redundante1 > 0)then
				if(Pn_QS_redundante1(1) > 0.0_dp)then
					t_LS_redundante1 = 1.0_dp /Pn_QS_redundante1(1)
				else
					t_LS_redundante1 = tMax
				endif
			endif
			
			do i1 = 1, nInfMax_redundante1
				if(Pn_QS_redundante1(i1) > 0.0_dp)then
					rho_medioQS_redundante1 = 1.0_dp * i1 * Pn_QS_redundante1(i1) + rho_medioQS_redundante1						! n medio
					rho2_medioQS_redundante1 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante1(i1) + rho2_medioQS_redundante1
					S_Schanon_redundante1 = S_Schanon_redundante1 - Pn_QS_redundante1(i1) * log(Pn_QS_redundante1(i1))
				endif
			enddo
			

			rho_medioQS_redundante1 = 1.0_dp * rho_medioQS_redundante1 / nInfMax_redundante1

			rho2_medioQS_redundante1 = 1.0_dp * rho2_medioQS_redundante1/ ((1.0_dp * nInfMax_redundante1)**2.0_dp)

			dev_rho_medioQS_redundante1 = (rho2_medioQS_redundante1 - rho_medioQS_redundante1**2.0_dp)**0.5_dp
			
			Xi_redundante1 = 1.0_dp * nInfMax_redundante1 * (rho2_medioQS_redundante1 - (rho_medioQS_redundante1**2.0_dp))/rho_medioQS_redundante1

!#######################################################################
!	Estatistica nos redundantes2
!#######################################################################
			
			Pn_QS_redundante2 = 1.0_dp * Pn_QS_redundante2 / sumPQS
			rho_medioQS_redundante2 = 0.0_dp
			rho2_medioQS_redundante2 = 0.0_dp
			S_Schanon_redundante2 = 0.0_dp		

			if(nInfMax_redundante2 > 0)then
				if(Pn_QS_redundante2(1) > 0.0_dp)then
					t_LS_redundante2 = 1.0_dp /Pn_QS_redundante2(1)
				else
					t_LS_redundante2 = tMax
				endif
			endif
			
			do i1 = 1, nInfMax_redundante2
				if(Pn_QS_redundante2(i1) > 0.0_dp)then
					rho_medioQS_redundante2 = 1.0_dp * i1 * Pn_QS_redundante2(i1) + rho_medioQS_redundante2						! n medio
					rho2_medioQS_redundante2 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante2(i1) + rho2_medioQS_redundante2
					S_Schanon_redundante2 = S_Schanon_redundante2 - Pn_QS_redundante2(i1) * log(Pn_QS_redundante2(i1))
				endif
			enddo
			

			rho_medioQS_redundante2 = 1.0_dp * rho_medioQS_redundante2 / nInfMax_redundante2

			rho2_medioQS_redundante2 = 1.0_dp * rho2_medioQS_redundante2/ ((1.0_dp * nInfMax_redundante2)**2.0_dp)

			dev_rho_medioQS_redundante2 = (rho2_medioQS_redundante2 - rho_medioQS_redundante2**2.0_dp)**0.5_dp
			
			Xi_redundante2 = 1.0_dp * nInfMax_redundante2 * (rho2_medioQS_redundante2 - (rho_medioQS_redundante2**2.0_dp))/rho_medioQS_redundante2
			
		end subroutine		

!#######################################################################
			
!#######################################################################		
		subroutine condicaoInicial_CG2(this, imuniza_tubos)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			class(grafo) :: this
			integer :: i1, j1, k1, lastCand, nInf_dummy
			logical, intent(in) :: imuniza_tubos
			

		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

				sigma = 0																			! Todos os nos sao suscetiveis
																				
				
				k1 = 0

lupy:			do i1 = 1, this%nodes

				if(lista_de_clusters(i1) /= i_comp_gigante)then
					sigma(i1) = 2
					cycle lupy
				endif	

				if(imuniza_tubos)then
					if(this%deg(i1) == 2)then
						sigma(i1) = 2
						cycle lupy
					endif
				endif
				
				k1 = k1 + 1
				infList(k1) = i1
				sigma(i1) = 1
				
			enddo lupy
						
			nInf = nInfMax

			nInf_n_redundante = nInfMax_n_redundante
			nInf_redundante1 = nInfMax_redundante1
			nInf_redundante2 = nInfMax_redundante2

			t = 0.0_dp;
			t_pos = 1;

		end subroutine

!#######################################################################

		subroutine condicaoInicial_CG(this)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			class(grafo) :: this
			integer :: i1, j1, k1, lastCand, nInf_dummy
			

		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

				sigma = 0																			! Todos os nos sao suscetiveis
																				
				
				k1 = 0

lupy:			do i1 = 1, this%nodes


!#######################################################################
!	Quem nao pertence aa CG, eu imunizo
				if(lista_de_clusters(i1) /= i_comp_gigante)then
					sigma(i1) = 2
					cycle lupy
				endif	
!#######################################################################
				
!#######################################################################
!	Quem pertence aa CG, comeca infectado.
				
				k1 = k1 + 1
				infList(k1) = i1
				sigma(i1) = 1
!#######################################################################

				
			enddo lupy

!#######################################################################
!	Eu setto o numero inifical de infectados, calculado no
!	Alloca_Medias...
						
			nInf = nInfMax

			nInf_n_redundante = nInfMax_n_redundante
			nInf_redundante = nInfMax_redundante
			
!#######################################################################

			
			t = 0._dp;
			t_pos = 1;

		end subroutine


!#######################################################################


		subroutine sisProcessRBS2(this,com_comp_gigante, seed)

!#######################################################################
			class(grafo), intent(in) :: this
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
			logical, intent(in) :: com_comp_gigante
!#######################################################################



!#######################################################################
!	nInfMax foi settado na rotina condicaoIncial_CG
!#######################################################################				 
				 if(lambda == 0.0_dp) lambda = dlambda
							
				 n_K = 0
				 
				 do i1 = 1, nInf
					n_K = n_K + this%deg(infList(i1))
				 enddo
				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp 
				Pn_QS_redundante1 = 0.0_dp 
				Pn_QS_redundante2 = 0.0_dp 
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		


			
loopDinamico:	do while(t <= tMax)
				
				
				
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = 1.0_dp * nInf * mu
				lambdaTotal = 1.0_dp * n_K * lambda
				rateTotal = 1.0_dp * muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				
	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt

	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################	

				
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante1(nInf_redundante1) = Pn_QS_redundante1(nInf_redundante1) + dt
					Pn_QS_redundante2(nInf_redundante2) = Pn_QS_redundante2(nInf_redundante2) + dt					
				endif

				
	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	

					if(nInf == 1)then
						sigma(infList(1)) = 0
						
						if(this%deg(infList(1)) > 2) nInf_n_redundante = 0				
						if(this%deg(infList(1)) == 1) nInf_redundante1 = 0
						if(this%deg(infList(1)) == 2) nInf_redundante2 = 0										

						n_K = 0							

						if(com_comp_gigante)then
							do
								i1 = gen%int(1,this%nodes)
								if(lista_de_clusters(i1) == i_comp_gigante) exit
							enddo
						else
							i1 = gen%int(1, this%nodes)
						endif
						
						infList(1) = i1	
						sigma(i1) = 1
						n_K = this%deg(i1)
						
						if(this%deg(i1) > 2) nInf_n_redundante = 1
						if(this%deg(i1) == 1) nInf_redundante1 = 1
						if(this%deg(i1) == 2) nInf_redundante2 = 1						
						cycle loopDinamico
					endif
			
	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
					
					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					
					if(this%deg(infList(i1)) > 2) nInf_n_redundante = nInf_n_redundante - 1
					if(this%deg(infList(i1)) == 1) nInf_redundante1 = nInf_redundante1 - 1
					if(this%deg(infList(i1)) == 2) nInf_redundante2 = nInf_redundante2 - 1																				
					
					infList(i1) = infList(nInf)
					nInf = nInf - 1

	!###################################################################################################################
	!	Como soh ha dois eventos possiveis, cura e infeccao, no caso de a cura nao ser sorteada, sobra a infeccao
	!	automaticamente.
	!	Um no infectado eh sorteado, estabelecemos a probabilidade dele infectar baseado no seu grau de infeccao
	!	probInf = 1.0_dp * lambda * this%deg(infList(i1))/ 
	!###################################################################################################################	


					
		   else	
loop_infecta:		do		

					i1 = gen%int(1, nInf)

					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle

						j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

						if(sigma(this%listAdj(j1)) == 0)then
							sigma(this%listAdj(j1)) = 1
							n_K = n_K + this%deg(this%listAdj(j1))
							nInf = nInf + 1																			! Numero de infectados aumenta.
							infList(nInf) = this%listAdj(j1)
											
							if(this%deg(this%listAdj(j1)) > 2) nInf_n_redundante = nInf_n_redundante + 1													
							if(this%deg(this%listAdj(j1)) == 1) nInf_redundante1 = nInf_redundante1 + 1
							if(this%deg(this%listAdj(j1)) == 2) nInf_redundante2 = nInf_redundante2 + 1
																											
						endif
						exit loop_infecta
			enddo loop_infecta
		   endif		
		   			
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para a Pn_QS, quanto para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################					
			
			write(*,*) "tamanho Pn_QS ", size(Pn_QS)
			
			Pn_QS = 1.0_dp * Pn_QS/sumPQS
									
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp			
			
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif

			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo

			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante

			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((1.0_dp * comp_gigante)**2.0_dp)

			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
						


!#######################################################################
!	Estatistica nos nao redundantes
!#######################################################################
			
			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp		
			
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS_n_redundante = 1.0_dp /Pn_QS_n_redundante(1)
			else
				t_LS_n_redundante = tMax
			endif

			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			

			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante

			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((1.0_dp * nInfMax_n_redundante)**2.0_dp)

			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante

!#######################################################################
!	Estatistica nos redundantes1
!#######################################################################
			
			Pn_QS_redundante1 = 1.0_dp * Pn_QS_redundante1 / sumPQS
			rho_medioQS_redundante1 = 0.0_dp
			rho2_medioQS_redundante1 = 0.0_dp
			S_Schanon_redundante1 = 0.0_dp		

			if(Pn_QS_redundante1(1) > 0.0_dp)then
				t_LS_redundante1 = 1.0_dp /Pn_QS_redundante1(1)
			else
				t_LS_redundante1 = tMax
			endif

			do i1 = 1, nInfMax_redundante1
				if(Pn_QS_redundante1(i1) > 0.0_dp)then
					rho_medioQS_redundante1 = 1.0_dp * i1 * Pn_QS_redundante1(i1) + rho_medioQS_redundante1						! n medio
					rho2_medioQS_redundante1 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante1(i1) + rho2_medioQS_redundante1
					S_Schanon_redundante1 = S_Schanon_redundante1 - Pn_QS_redundante1(i1) * log(Pn_QS_redundante1(i1))
				endif
			enddo
			

			rho_medioQS_redundante1 = 1.0_dp * rho_medioQS_redundante1 / nInfMax_redundante1

			rho2_medioQS_redundante1 = 1.0_dp * rho2_medioQS_redundante1/ ((1.0_dp * nInfMax_redundante1)**2.0_dp)

			dev_rho_medioQS_redundante1 = (rho2_medioQS_redundante1 - rho_medioQS_redundante1**2.0_dp)**0.5_dp
			
			Xi_redundante1 = 1.0_dp * nInfMax_redundante1 * (rho2_medioQS_redundante1 - (rho_medioQS_redundante1**2.0_dp))/rho_medioQS_redundante1

!#######################################################################
!	Estatistica nos redundantes2
!#######################################################################
			
			Pn_QS_redundante2 = 1.0_dp * Pn_QS_redundante2 / sumPQS
			rho_medioQS_redundante2 = 0.0_dp
			rho2_medioQS_redundante2 = 0.0_dp
			S_Schanon_redundante2 = 0.0_dp		

			if(Pn_QS_redundante2(1) > 0.0_dp)then
				t_LS_redundante2 = 1.0_dp /Pn_QS_redundante2(1)
			else
				t_LS_redundante2 = tMax
			endif

			do i1 = 1, nInfMax_redundante2
				if(Pn_QS_redundante2(i1) > 0.0_dp)then
					rho_medioQS_redundante2 = 1.0_dp * i1 * Pn_QS_redundante2(i1) + rho_medioQS_redundante2						! n medio
					rho2_medioQS_redundante2 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante2(i1) + rho2_medioQS_redundante2
					S_Schanon_redundante2 = S_Schanon_redundante2 - Pn_QS_redundante2(i1) * log(Pn_QS_redundante2(i1))
				endif
			enddo
			

			rho_medioQS_redundante2 = 1.0_dp * rho_medioQS_redundante2 / nInfMax_redundante2

			rho2_medioQS_redundante2 = 1.0_dp * rho2_medioQS_redundante2/ ((1.0_dp * nInfMax_redundante2)**2.0_dp)

			dev_rho_medioQS_redundante2 = (rho2_medioQS_redundante2 - rho_medioQS_redundante2**2.0_dp)**0.5_dp
			
			Xi_redundante2 = 1.0_dp * nInfMax_redundante2 * (rho2_medioQS_redundante2 - (rho_medioQS_redundante2**2.0_dp))/rho_medioQS_redundante2


			
		end subroutine		

!#######################################################################
!	Esse eh o mais atual

		subroutine sisProcessRBS1(this,com_comp_gigante, seed)

!#######################################################################
			class(grafo), intent(in) :: this
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
			logical, intent(in) :: com_comp_gigante
!#######################################################################



!#######################################################################
!	nInfMax foi settado na rotina condicaoIncial_CG
!#######################################################################				 
				 if(lambda == 0.0_dp) lambda = dlambda
							
				 n_K = 0
				 
				 do i1 = 1, nInf
					n_K = n_K + this%deg(infList(i1))
				 enddo
				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp 
				Pn_QS_redundante1 = 0.0_dp 
				Pn_QS_redundante2 = 0.0_dp 
				
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		


			
loopDinamico:	do while(t <= tMax)
				
				
				
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = 1.0_dp * nInf * mu
				lambdaTotal = 1.0_dp * n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				
	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt

	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################	
				
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante1(nInf_redundante1) = Pn_QS_redundante1(nInf_redundante1) + dt
					Pn_QS_redundante2(nInf_redundante2) = Pn_QS_redundante2(nInf_redundante2) + dt					
				endif


				
	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	

					if(nInf == 1)then
						sigma(infList(1)) = 0
						
						if(this%deg(infList(1)) > 2) nInf_n_redundante = 0				
						if(this%deg(infList(1)) == 1) nInf_redundante1 = 0
						if(this%deg(infList(1)) == 2) nInf_redundante2 = 0										

						n_K = 0							

						if(com_comp_gigante)then
							do
								i1 = gen%int(1,this%nodes)
								if(lista_de_clusters(i1) == i_comp_gigante) exit
							enddo
						else
							i1 = gen%int(1, this%nodes)
						endif
						
						infList(1) = i1	
						sigma(i1) = 1
						n_K = this%deg(i1)
						
						if(this%deg(i1) > 2) nInf_n_redundante = 1
						if(this%deg(i1) == 1) nInf_redundante1 = 1
						if(this%deg(i1) == 2) nInf_redundante2 = 1
												
						cycle loopDinamico
					endif
			
	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
					
					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					
					if(this%deg(infList(i1)) > 2) nInf_n_redundante = nInf_n_redundante - 1
					if(this%deg(infList(i1)) == 1) nInf_redundante1 = nInf_redundante1 - 1
					if(this%deg(infList(i1)) == 2) nInf_redundante2 = nInf_redundante2 - 1															
					
					infList(i1) = infList(nInf)
					nInf = nInf - 1

	!###################################################################################################################
	!	Como soh ha dois eventos possiveis, cura e infeccao, no caso de a cura nao ser sorteada, sobra a infeccao
	!	automaticamente.
	!	Um no infectado eh sorteado, estabelecemos a probabilidade dele infectar baseado no seu grau de infeccao
	!	probInf = 1.0_dp * lambda * this%deg(infList(i1))/ 
	!###################################################################################################################	


					
		   else	
loop_infecta:		do		

					i1 = gen%int(1, nInf)

					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle

						j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

						if(sigma(this%listAdj(j1)) == 0)then
							sigma(this%listAdj(j1)) = 1
							n_K = n_K + this%deg(this%listAdj(j1))
							
							nInf = nInf + 1																			! Numero de infectados aumenta.
							infList(nInf) = this%listAdj(j1)
											
							if(this%deg(this%listAdj(j1)) > 2) nInf_n_redundante = nInf_n_redundante + 1													
							if(this%deg(this%listAdj(j1)) == 1) nInf_redundante1 = nInf_redundante1 + 1
							if(this%deg(this%listAdj(j1)) == 2) nInf_redundante2 = nInf_redundante2 + 1
																											
						endif
						exit loop_infecta
			enddo loop_infecta
		   endif		
		   			
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para a Pn_QS, quanto para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################					
							
			do i1 = 1, size(Pn_QS)			
				Pn_QS(i1) = 1.0_dp * Pn_QS(i1) / sumPQS
			enddo
									
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp			
			
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif

			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo

			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante

			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((1.0_dp * comp_gigante)**2.0_dp)

			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
						


!#######################################################################
!	Estatistica nos nao redundantes
!#######################################################################
			
			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp		
			
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS_n_redundante = 1.0_dp /Pn_QS_n_redundante(1)
			else
				t_LS_n_redundante = tMax
			endif

			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			

			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante

			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((1.0_dp * nInfMax_n_redundante)**2.0_dp)

			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante

!#######################################################################
!	Estatistica nos redundantes1
!#######################################################################
			
			Pn_QS_redundante1 = 1.0_dp * Pn_QS_redundante1 / sumPQS
			rho_medioQS_redundante1 = 0.0_dp
			rho2_medioQS_redundante1 = 0.0_dp
			S_Schanon_redundante1 = 0.0_dp		

			if(nInfMax_redundante1 > 0)then
				if(Pn_QS_redundante1(1) > 0.0_dp)then
					t_LS_redundante1 = 1.0_dp /Pn_QS_redundante1(1)
				else
					t_LS_redundante1 = tMax
				endif
			endif
			
			do i1 = 1, nInfMax_redundante1
				if(Pn_QS_redundante1(i1) > 0.0_dp)then
					rho_medioQS_redundante1 = 1.0_dp * i1 * Pn_QS_redundante1(i1) + rho_medioQS_redundante1						! n medio
					rho2_medioQS_redundante1 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante1(i1) + rho2_medioQS_redundante1
					S_Schanon_redundante1 = S_Schanon_redundante1 - Pn_QS_redundante1(i1) * log(Pn_QS_redundante1(i1))
				endif
			enddo
			

			rho_medioQS_redundante1 = 1.0_dp * rho_medioQS_redundante1 / nInfMax_redundante1

			rho2_medioQS_redundante1 = 1.0_dp * rho2_medioQS_redundante1/ ((1.0_dp * nInfMax_redundante1)**2.0_dp)

			dev_rho_medioQS_redundante1 = (rho2_medioQS_redundante1 - rho_medioQS_redundante1**2.0_dp)**0.5_dp
			
			Xi_redundante1 = 1.0_dp * nInfMax_redundante1 * (rho2_medioQS_redundante1 - (rho_medioQS_redundante1**2.0_dp))/rho_medioQS_redundante1

!#######################################################################
!	Estatistica nos redundantes2
!#######################################################################
			
			Pn_QS_redundante2 = 1.0_dp * Pn_QS_redundante2 / sumPQS
			rho_medioQS_redundante2 = 0.0_dp
			rho2_medioQS_redundante2 = 0.0_dp
			S_Schanon_redundante2 = 0.0_dp		

			if(nInfMax_redundante2 > 0)then
				if(Pn_QS_redundante2(1) > 0.0_dp)then
					t_LS_redundante2 = 1.0_dp /Pn_QS_redundante2(1)
				else
					t_LS_redundante2 = tMax
				endif
			endif
			
			do i1 = 1, nInfMax_redundante2
				if(Pn_QS_redundante2(i1) > 0.0_dp)then
					rho_medioQS_redundante2 = 1.0_dp * i1 * Pn_QS_redundante2(i1) + rho_medioQS_redundante2						! n medio
					rho2_medioQS_redundante2 =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante2(i1) + rho2_medioQS_redundante2
					S_Schanon_redundante2 = S_Schanon_redundante2 - Pn_QS_redundante2(i1) * log(Pn_QS_redundante2(i1))
				endif
			enddo
			

			rho_medioQS_redundante2 = 1.0_dp * rho_medioQS_redundante2 / nInfMax_redundante2

			rho2_medioQS_redundante2 = 1.0_dp * rho2_medioQS_redundante2/ ((1.0_dp * nInfMax_redundante2)**2.0_dp)

			dev_rho_medioQS_redundante2 = (rho2_medioQS_redundante2 - rho_medioQS_redundante2**2.0_dp)**0.5_dp
			
			Xi_redundante2 = 1.0_dp * nInfMax_redundante2 * (rho2_medioQS_redundante2 - (rho_medioQS_redundante2**2.0_dp))/rho_medioQS_redundante2
			
		end subroutine		


!####################################################################################
!		Processo SIS com condicao de contorno reflexiva
!####################################################################################
		subroutine sisProcessRBS(this,com_comp_gigante, seed)

!#######################################################################
			class(grafo), intent(in) :: this
			integer :: seed
			integer :: i1, j1, k1
			real(dp) :: probInf, prob_dist_t
			real(dp) :: propRej, startCura, finishCura, timeCura, timeCuraAc, propCura, startLat,&
             finishLat, timeLat, timeLatAc, propLat
			logical, intent(in) :: com_comp_gigante
!#######################################################################



!#######################################################################
!	nInfMax foi settado na rotina condicaoIncial_CG
!#######################################################################				 
				 if(lambda == 0.0_dp) lambda = dlambda
							
				 n_K = 0
				 
				 do i1 = 1, nInf
					n_K = n_K + this%deg(infList(i1))
				 enddo
				 
	!###################################################################################################################
	!			Aqui comeca o regime quasi-estacionario, passada a fase transiente ou de relaxacao.
	!###################################################################################################################	

				Pn_QS = 0.0_dp
				Pn_QS_n_redundante = 0.0_dp 
				Pn_QS_redundante = 0.0_dp 
	!###################################################################
	!	Abre os arquivos para salvar a dinamica
	!###################################################################		


			
loopDinamico:	do while(t <= tMax)
				
				
				
	!###################################################################################################################
	!	Estabelecemos todas as taxas de eventos aqui. Soh pode se curar quem estah infectado, portanto, a taxa
	!	de cura eh muTotal = nInf * mu. Como a infeccao soh pode fluir atraves das conexoes com nos infectados,
	!	a taxa de infeccao total eh lambdaTotal = n_K * lambda. Notar que ela flui mesmo se o vizinho for infectado,
	!	resultando aih num processo fantasma.
	!	A taxa total de eventos eh a soma dos dois anteriores
	!###################################################################################################################	
				
				muTotal = nInf * mu
				lambdaTotal = n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				
	!###################################################################################################################
	!	O intervalo entre dois eventos eh calculado por meio de uma distribuicao de probabilidades em decaimento
	!	exponencial e o tempo eh atualizado.
	!###################################################################################################################	
				prob_dist_t = gen%rnd()
				dt = -1.0_dp * log(max(1e-12,prob_dist_t))/rateTotal
				t = t + dt

	!###################################################################################################################
	!	Passado o tempo de relaxacao, passamos a calcular a probabilidade quasi-estacionaria.
	!###################################################################################################################	

				
				if(t >= t_relax)then
					Pn_QS(nInf) = Pn_QS(nInf) + dt
					Pn_QS_n_redundante(nInf_n_redundante) = Pn_QS_n_redundante(nInf_n_redundante) + dt
					Pn_QS_redundante(nInf_redundante) = Pn_QS_redundante(nInf_redundante) + dt
				endif

				
	!###############################################################################################################
	!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
	!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
	!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
	!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
	!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
	!###############################################################################################################

				
		   prob = gen%rnd()
				
		   if(prob <= m)then


	!####################################################################################################################
	!	Aqui ocorre a sutileza da Condição de Contorno Reflexiva. Caso o numero de nos infectados seja 1,
	!	nao permitimos que haja cura. No entanto, para evitar a perda de ergodicidade, curamos o antigo no infectado, 
	!	retiramos seus stubs infectados do pool de stubs infectados, que passa a ser 0 agora, sorteamos um no,
	!	colocamos o no no local do antigo no infectado na lista de infectados, retirando assim o antigo da lista,
	!	e infectamos o no sorteado, mudando seu estado na lista sigma. Entao adicionamos os stubs agora infectados
	!	do novo no ao pool de stubs infectados. Entao voltamos a correr a dinamica.
	!####################################################################################################################	

					if(nInf == 1)then
						sigma(infList(1)) = 0
						
						if(this%deg(infList(1)) >= 3) nInf_n_redundante = 0				
						if(this%deg(infList(1)) < 3) nInf_redundante = 0				

						n_K = 0							

						if(com_comp_gigante)then
							do
								i1 = gen%int(1,this%nodes)
								if(lista_de_clusters(i1) == i_comp_gigante) exit
							enddo
						else
							i1 = gen%int(1, this%nodes)
						endif
						
						infList(1) = i1	
						sigma(i1) = 1
						n_K = this%deg(i1)
						
						if(this%deg(i1) >= 3) nInf_n_redundante = 1
						if(this%deg(i1) < 3) nInf_redundante = 1						
						cycle loopDinamico
					endif
			
	!###################################################################################################################
	!	Caso possa haver cura, sorteamos um no infectado, mudamos seu estado para curado, na lista sigma,
	!	subtraimos seus stubs do pool de stubs infectados, colocamos o ultimo no da lista no seu local e diminuimos
	!	o numero de nos infectados em uma unidade.
	!###################################################################################################################	
					
					i1 = gen%int(1, nInf)
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					
					if(this%deg(infList(i1)) >= 3) nInf_n_redundante = nInf_n_redundante - 1
					if(this%deg(infList(i1)) < 3) nInf_redundante = nInf_redundante - 1															
					
					infList(i1) = infList(nInf)
					nInf = nInf - 1

	!###################################################################################################################
	!	Como soh ha dois eventos possiveis, cura e infeccao, no caso de a cura nao ser sorteada, sobra a infeccao
	!	automaticamente.
	!	Um no infectado eh sorteado, estabelecemos a probabilidade dele infectar baseado no seu grau de infeccao
	!	probInf = 1.0_dp * lambda * this%deg(infList(i1))/ 
	!###################################################################################################################	


					
		   else	
loop_infecta:		do		

					i1 = gen%int(1, nInf)

					probInf = 1.0_dp * this%deg(infList(i1))/this%degMax
					prob = gen%rnd()
					
					if(prob	> probInf)cycle

						j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)																													! de um sitio

						if(sigma(this%listAdj(j1)) == 0)then
							sigma(this%listAdj(j1)) = 1
							n_K = n_K + this%deg(this%listAdj(j1))
							nInf = nInf + 1																			! Numero de infectados aumenta.
							infList(nInf) = this%listAdj(j1)
											
							if(this%deg(this%listAdj(j1)) >= 3) nInf_n_redundante = nInf_n_redundante + 1													
							if(this%deg(this%listAdj(j1)) < 3) nInf_redundante = nInf_redundante + 1																				
						endif
						exit loop_infecta
			enddo loop_infecta
		   endif		
		   			
		enddo loopDinamico


!#############################################################################################################################
!	Aqui nos atualizamos todas as variaveis dinamicas
!#############################################################################################################################
			
			!###########################################################
			!	Serve tanto para a Pn_QS, quanto para Pn_QS_n_redundante	
				sumPQS = sum(Pn_QS)
			!###########################################################					
				
						
			Pn_QS = 1.0_dp * Pn_QS / sumPQS
									
			rho_medioQS = 0.0_dp
			rho2_medioQS = 0.0_dp
			S_Schanon = 0.0_dp			
			
			if(Pn_QS(1) > 0.0_dp)then
				t_LS = 1.0_dp /Pn_QS(1)
			else
				t_LS = tMax
			endif

			
			do i1 = 1, nInfMax
				if(Pn_QS(i1) > 0.0_dp)then
					rho_medioQS = 1.0_dp * i1 * Pn_QS(i1) + rho_medioQS						! n medio
					rho2_medioQS =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS(i1) + rho2_medioQS
					S_Schanon = S_Schanon - Pn_QS(i1) * log(Pn_QS(i1))
				endif
			enddo

			rho_medioQS = 1.0_dp * rho_medioQS / comp_gigante

			rho2_medioQS = 1.0_dp * rho2_medioQS/ ((comp_gigante)**2.0_dp)

			dev_rho_medioQS = (rho2_medioQS - rho_medioQS**2.0_dp)**0.5_dp
			
			Xi = 1.0_dp * comp_gigante * (rho2_medioQS - (rho_medioQS**2.0_dp))/rho_medioQS
						


!#######################################################################
!	Estatistica nos nao redundantes
!#######################################################################
			
			Pn_QS_n_redundante = 1.0_dp * Pn_QS_n_redundante / sumPQS
			rho_medioQS_n_redundante = 0.0_dp
			rho2_medioQS_n_redundante = 0.0_dp
			S_Schanon_n_redundante = 0.0_dp		
			
			if(Pn_QS_n_redundante(1) > 0.0_dp)then
				t_LS_n_redundante = 1.0_dp /Pn_QS_n_redundante(1)
			else
				t_LS_n_redundante = tMax
			endif

			do i1 = 1, nInfMax_n_redundante
				if(Pn_QS_n_redundante(i1) > 0.0_dp)then
					rho_medioQS_n_redundante = 1.0_dp * i1 * Pn_QS_n_redundante(i1) + rho_medioQS_n_redundante						! n medio
					rho2_medioQS_n_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_n_redundante(i1) + rho2_medioQS_n_redundante
					S_Schanon_n_redundante = S_Schanon_n_redundante - Pn_QS_n_redundante(i1) * log(Pn_QS_n_redundante(i1))
				endif
			enddo
			

			rho_medioQS_n_redundante = 1.0_dp * rho_medioQS_n_redundante / nInfMax_n_redundante

			rho2_medioQS_n_redundante = 1.0_dp * rho2_medioQS_n_redundante/ ((nInfMax_n_redundante)**2.0_dp)

			dev_rho_medioQS_n_redundante = (rho2_medioQS_n_redundante - rho_medioQS_n_redundante**2.0_dp)**0.5_dp
			
			Xi_n_redundante = 1.0_dp * nInfMax_n_redundante * (rho2_medioQS_n_redundante - (rho_medioQS_n_redundante**2.0_dp))/rho_medioQS_n_redundante

!#######################################################################
!	Estatistica nos redundantes
!#######################################################################
			
			Pn_QS_redundante = 1.0_dp * Pn_QS_redundante / sumPQS
			rho_medioQS_redundante = 0.0_dp
			rho2_medioQS_redundante = 0.0_dp
			S_Schanon_redundante = 0.0_dp		

			if(nInfMax_redundante > 0)then
				if(Pn_QS_redundante(1) > 0.0_dp)then
					t_LS_redundante = 1.0_dp /Pn_QS_redundante(1)
				else
					t_LS_redundante = tMax
				endif
			endif
			
			do i1 = 1, nInfMax_redundante
				if(Pn_QS_redundante(i1) > 0.0_dp)then
					rho_medioQS_redundante = 1.0_dp * i1 * Pn_QS_redundante(i1) + rho_medioQS_redundante						! n medio
					rho2_medioQS_redundante =  ((1.0_dp *i1) ** 2.0_dp) * Pn_QS_redundante(i1) + rho2_medioQS_redundante
					S_Schanon_redundante = S_Schanon_redundante - Pn_QS_redundante(i1) * log(Pn_QS_redundante(i1))
				endif
			enddo
			

			rho_medioQS_redundante = 1.0_dp * rho_medioQS_redundante / nInfMax_redundante

			rho2_medioQS_redundante = 1.0_dp * rho2_medioQS_redundante/ ((nInfMax_redundante)**2.0_dp)

			dev_rho_medioQS_redundante = (rho2_medioQS_redundante - rho_medioQS_redundante**2.0_dp)**0.5_dp
			
			Xi_redundante = 1.0_dp * nInfMax_redundante * (rho2_medioQS_redundante - (rho_medioQS_redundante**2.0_dp))/rho_medioQS_redundante

			
		end subroutine		
!#######################################################################
	
!#######################################################################		
 
 		
		subroutine condicaoInicial(this)
		
			!###############################################################################
			!	Variaveis auxiliares
			!###############################################################################
			
			integer :: i1, j1, k1, lastCand
			
			
			!###############################################################################
			!	Variaveis inicializacao da dinamica
			!###############################################################################
			
			
			
			!###############################################################################
			!	A rede, ou substrato
			!###############################################################################

			class(grafo), intent(in) :: this

			!###############################################################################
			!	Lista aulixiar
			!###############################################################################
	

			
		!#######################################################################################
		!			Situacao inicial da dinamica
		!#######################################################################################

			sigma = 0																			! Todos os nos sao suscetiveis
			nInf = int(nInf0 * comp_gigante)		 												! Quantidade de nos inicialmente infectados
			nSus = comp_gigante - nInf															! Quantidade de nos inicialmente suscetiveis
											
!			write(*,*) "Numero inicial de infectados ", nInf
										
			if(nInf0 < 1._dp)then
				if(allocated(listAux)) deallocate(listAux)	
				allocate(listAux(nInf))
			endif
							
			if(nInf0 < 1._dp)then
					j1 = 0																																		
					do i1 = 1, this%nodes
						if(lista_de_clusters(i1) == i_comp_gigante)then
							j1 = j1 + 1
							listAux(j1) = i1
							if(j1 == nInf) exit
						endif
					enddo
					k1 = 1;	lastCand = this%nodes
					do while(k1 <= nInf)

						i1 = gen%int(1, lastCand)												! lastCand pertence a geraRede
						infList(k1) = listAux(i1)
						sigma(listAux(i1)) = 1

						listAux(i1) = listAux(lastCand)
						lastCand = lastCand - 1

						k1 = k1 + 1
					enddo
			else
					do i1 = 1, this%nodes
						infList(i1) = i1
						sigma(i1) = 1
					enddo
			endif											

			t = 0._dp;
			t_pos = 1;
		end subroutine
		
!####################################################################### 
 
 
		subroutine sisProcess(this)


			!###################################################################################
			!	A rede, ou substrato
			!###################################################################################

			class(grafo), intent(in) :: this
			integer ::  n_K
			integer :: i1, j1, k1
			real(dp) :: probInf
			 
			 n_K = 0
			 
			 do i1 = 1, nInf
				n_K = n_K + this%deg(infList(i1))
			 enddo
			 
			 
			
loopDinamico:	do while(t <= tMax)
				
				muTotal = nInf * mu
				lambdaTotal = n_K * lambda
				rateTotal = muTotal + lambdaTotal
				m = 1.0_dp * muTotal/rateTotal
				
				dt = -1.0_dp * log(1._dp - gen%rnd())/rateTotal
				t = t + dt
				
				prob = gen%rnd()
				if(prob <= m)then																	! Se o processo da vez eh cura.
					i1 = gen%int(1, nInf)
					
					sigma(infList(i1)) = 0
					n_K = n_K - this%deg(infList(i1))
					infList(i1) = infList(nInf)
					nInf = nInf - 1		
!					write(*,*) "Cura! N Infectados = ", nInf
				else	
loop_infecta:				do														! Se o processo da vez eh infeccao.
						i1 = gen%int(1, nInf)
						probInf = 1._dp * this%deg(infList(i1))/this%degMax
						prob = gen%rnd()
						if(prob	<= probInf)then
							j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)
																		! Posicao na lista de adjacencia correspondente aos vizinhos
																											! de um sitio
							if(sigma(this%listAdj(j1)) == 0)then								
								sigma(this%listAdj(j1)) = 1
								n_K = n_K + this%deg(this%listAdj(j1))													! Muda o estado do sitio.
								nInf = nInf + 1													! Numero de infectados aumenta.
								infList(nInf) = this%listAdj(j1)												! Coloco esse novo infectado na primeira posicao disponivel
																											! na lista de infectados.
							endif
		!					write(*,*) "Infeccao! Numero de Infectados= ", nInf
							exit loop_infecta
						endif
					enddo loop_infecta
				endif									
				
				do while(t_pos <= t)
					rho_medio(t_pos) = rho_medio(t_pos) + 1._dp * nInf/this%nodes
					t_medio(t_pos) = t_medio(t_pos) + t
					t_samp(t_pos) = t_samp(t_pos) + 1
					
					if(nInf .ne. 0)then
						t_MaxVis = max(t_MaxVis, t_pos)												! Qual o tempo recorde que uma infeccao
					endif																			! sobreviveu ate agora na rede?	
					t_pos = t_pos + 1
				enddo

				if(nInf == 0) exit loopDinamico														! Porque o estado absorvente foi alcancado

			enddo loopDinamico			
		end subroutine
		
!#######################################################################

		
!####################################################################################
!		Processo de Contato
!####################################################################################

	
			subroutine contactProcess(this)


			!###################################################################################
			!	A rede, ou substrato
			!###################################################################################

			class(grafo), intent(in) :: this
			integer :: i1, j1, k1
			
			
			m = 1.0_dp * mu/(lambda + mu)
			
loopDinamico:	do while(t <= tMax)
				
				muTotal = nInf * mu
				lambdaTotal = nInf * lambda
				rateTotal = muTotal + lambdaTotal
				
				dt = -1.0_dp * log(1._dp - gen%rnd())/rateTotal
				t = t + dt
				
		!###############################################################################################################
		!	Testa as chances de ocorrer cura. Se sim, um no eh escolhido dentre os nos infectados,
		!	obtidos da lista de infectados, seu estado eh modificado de 1 para 0, na lista de estados sigma,
		!	ele eh retirado da lista de infectados, atraves da adicao do ultimo noh da lista em seu local,
		!	e o tamanho virtual da lista eh subtraido de uma unidade. Se ele proprio for o ultimo da lista,
		!	a subtracao de uma unidade do tamanho da lista jah o retira da lista
		!###############################################################################################################
				prob = gen%rnd()
				if(prob <= m)then																							! Se o processo da vez eh cura.
					i1 = gen%int(1, nInf)
					
					sigma(infList(i1)) = 0
					infList(i1) = infList(nInf)
					nInf = nInf - 1		
!					write(*,*) "Cura! N Infectados = ", nInf
				else																				! Se o processo da vez eh infeccao.
					
					i1 = gen%int(1, nInf)	
					j1 = gen%int(this%aux(infList(i1)), this%aux(infList(i1)) + this%deg(infList(i1)) - 1)							! Posicao na lista de adjacencia correspondente aos vizinhos
																									! de um sitio
					if(sigma(this%listAdj(j1)) == 0)then								
						sigma(this%listAdj(j1)) = 1													! Muda o estado do sitio.
						nInf = nInf + 1																! Numero de infectados aumenta.
						infList(nInf) = this%listAdj(j1)												! Coloco esse novo infectado na primeira posicao disponivel
																									! na lista de infectados.
					endif
!					write(*,*) "Infeccao! Numero de Infectados= ", nInf
				endif									
				
				do while(t_pos <= t)
					rho_medio(t_pos) = rho_medio(t_pos) + 1.0_dp * nInf/this%nodes
					t_medio(t_pos) = t_medio(t_pos) + t
					t_samp(t_pos) = t_samp(t_pos) + 1
					
					if(nInf .ne. 0)then
						t_MaxVis = max(t_MaxVis, t_pos)												! Qual o tempo recorde que uma infeccao
					endif																			! sobreviveu ate agora na rede?	
					t_pos = t_pos + 1
				enddo

				if(nInf == 0) exit loopDinamico														! Porque o estado absorvente foi alcancado

			enddo loopDinamico			
		end subroutine
						
		
		
!####################################################################################
!			Esta subrotina calcula as medias da epidemia
!####################################################################################
		
		subroutine calculaMedias(this,seed)
			integer :: seed
			integer :: j1, j2
			class(grafo), intent(in) :: this


			!###################################################
			!	Fase de inicializacao da subrotina
			!###################################################
			
			call gen%init(seed)
			
			
			
			do j1 = 1, nSamples
				write(*,*) "Calculando amostra de numero ", j1
				call gastaRnd()
				call condicaoInicial(this)
				call sisProcess(this)
				
			enddo
			
			open(unit=10, file='t_vs_rho.dat', status='unknown')
			do j1 = 1, t_MaxVis
					write(10, *) 1.0_dp * t_medio(j1)/t_samp(j1), ",", 1.0_dp * rho_medio(j1)/nSamples
			enddo
			close(10)

			
			!#############################################
			!	Subrotinas internas
			!#############################################
			
			contains
				subroutine gastaRnd()
					integer :: j4
					real(dp) :: j5
					do j4 = 1, 10 * nSamples
						j5 = gen%rnd()
					enddo
				end subroutine	
				
		end subroutine
		


!#########################################################################################
!	Classifica o tamanho da componenten gigante dinamica
!#########################################################################################

!####################################################################################	

	subroutine sub_classifica_clusters_nao_imune(this, criaArquivo, label, arquivo)

		class(grafo) :: this
		integer, intent(in) :: label
		character(len = *) :: arquivo
		logical :: criaArquivo
		real(dp) ::  grauMedio2
		integer :: sumgrau2
	!############################################################
	!	Variaveis mudas
	!############################################################
		integer :: i1, i2, i3, i4
	!############################################################
	!	Listas auxiliares internas
	!############################################################
		
		integer, allocatable :: lista_tam_clusters_din(:), vertices_a_testar(:), rede_original(:)

	!############################################################
	!	Variaveis internas
	!############################################################
		integer :: n_zeros, indice_cluster, novidade, agregados, pos_novosVertices, pos_ultimosVertices




!#################################################################################################################################################	
!		Alocamento de listas auxiliares internas
!#################################################################################################################################################
	
	write(*,*) "Rotina de clusteamento dinamico chamada com sucesso"
	
	
	if(allocated(rede_original)) deallocate(rede_original)
		allocate(rede_original(this%nodes))

	do i1=1, this%nodes
		rede_original(i1) = i1
	enddo


	if(allocated(vertices_a_testar)) deallocate(vertices_a_testar)
		allocate(vertices_a_testar(this%nodes))

		vertices_a_testar = 0


!#################################################################################################################################################	
!		Alocamento de listas importantes
!#################################################################################################################################################


	if(allocated(lista_de_clusters_din)) deallocate(lista_de_clusters_din)
		allocate(lista_de_clusters_din(this%nodes))
		lista_de_clusters_din = 0
		
		
	if(allocated(lista_tam_clusters_din)) deallocate(lista_tam_clusters_din)
		allocate(lista_tam_clusters_din(this%nodes))

		lista_tam_clusters_din = 0

		

	!###################################################################################################
	!				Inicializacoes importantes
	!########################################this%nodes###########################################################

!#######################################################################
			n_zeros = 0
			indice_cluster = 1			!O primeiro cluster tem indice 1.						

ch_first:	do i1 = 1, this%nodes
				if(rede_original(i1) /= 0)then
					
					n_zeros = n_zeros + 1
					rede_original(i1) = 0
					
					if(sigma(i1) == 2) cycle
					
					vertices_a_testar(1) = i1
					lista_tam_clusters_din(indice_cluster) = lista_tam_clusters_din(indice_cluster) + 1
					lista_de_clusters_din(vertices_a_testar(1)) = indice_cluster									
					exit ch_first	
					
				endif
			enddo ch_first
!#######################################################################
			
			
				
			pos_novosVertices = 1
			pos_ultimosVertices = 1


		
		do while(n_zeros < this%nodes)
		
			novidade = 1				! Necessario inicializar esse cara como nao nulo
		

			do while(novidade /= 0)
						
															 				
					agregados = pos_ultimosVertices
				
					novidade = 0
				
					do i1 = pos_novosVertices, pos_ultimosVertices
						
procura_agregados:				do i2 = this%aux(vertices_a_testar(i1)), this%aux(vertices_a_testar(i1)) + this%deg(vertices_a_testar(i1)) - 1
												
									if(rede_original(this%listAdj(i2)) /= 0)then

										rede_original(this%listAdj(i2)) = 0									
										n_zeros = n_zeros + 1
										
										if(sigma(this%listAdj(i2)) == 2) cycle procura_agregados
										
										agregados = agregados + 1
								
										vertices_a_testar(agregados) = this%listAdj(i2)
													
										lista_de_clusters_din(this%listAdj(i2)) = indice_cluster
						
										lista_tam_clusters_din(indice_cluster) = lista_tam_clusters_din(indice_cluster) + 1					
											
										novidade = 1		
															!Aqui diz 'Hey! temos alguem
															! para nos apresentar novos vizinhos!
															! como novidade' no cluster.						
									endif
								enddo procura_agregados
					enddo

						pos_novosVertices = pos_ultimosVertices + 1
						pos_ultimosVertices = agregados
			enddo
		
			pos_ultimosVertices = pos_novosVertices		! Quando nada de novo acontece, eh porque fechou um cluster
														! Daih, passamos a olhar  adiante na lista de vertices a testar
														! O numero de agregados eh o numero que representa
														! O tamanho do cluster anterior
		loop_p:	do i2 = 1, this%nodes
				
				if(rede_original(i2) /= 0) then

					rede_original(i2) = 0									! Antes eu zerava esse cara			
					n_zeros = n_zeros + 1

					if(sigma(i2) == 2) cycle loop_p
					                     
					indice_cluster = indice_cluster + 1						!Quer dizer que ha um novo cluster
					
					vertices_a_testar(pos_novosVertices) = i2				! Como nao zero mais, o novo a testar
																			 
					lista_de_clusters_din(i2) = indice_cluster
																							!					
					lista_tam_clusters_din(indice_cluster) = lista_tam_clusters_din(indice_cluster) + 1
								
					exit  loop_p
			
				endif
			enddo loop_p
		enddo

!#######################################################################
!	Escreve num arquivo ou nao?
!#######################################################################
	
		if(criaArquivo) then
			call abreArquivo(label, "sitio_e_tamanho_cluster.dat")

			do i1 =1, this%nodes
				if(lista_tam_clusters_din(i1) /= 0) then
					write(label,*) i1, lista_tam_clusters_din(i1)
				endif
			enddo

			close(label)

			call abreArquivo(label + 1, arquivo)
			write(label + 1, *) "Vertice ", "GrauDin", "Cluster"
		
			do i1 = 1, size(lista_de_clusters_din)
				write(label + 1,*) i1, this%deg(i1), lista_de_clusters_din(i1)
			enddo
		
			close(label + 1)
		endif

!#######################################################################
!	Faz a estatistica dos clusters
!#######################################################################
	
		comp_gigante_din = maxval(lista_tam_clusters_din)

		write(*,*) "ao que parece, comp_gigante_din e : ", comp_gigante_din
		
		do i1 = 1, size(lista_tam_clusters_din)
			if(lista_tam_clusters_din(i1) == comp_gigante_din)then
				
				i_comp_gigante_din = i1
				exit
				
			endif
		enddo
	
		write(*,*) "ao que parece, i_comp_gigante_din e : ", i_comp_gigante_din

!#######################################################################
!	Vamos testar se ha apenas uma comp gigante
!#######################################################################


		i2 = 0
		do i1 = 1, size(lista_tam_clusters_din)
			if(lista_tam_clusters_din(i1) == comp_gigante_din) i2 = i2 + 1 
		enddo
		
		if( i2 > 1 ) stop "Ha mais de uma componente gigante"


!#######################################################################
!	Testa se o grau maximo continua o mesmo
!#######################################################################

		write(*,*) 'O grau maximo era: ', maxval(this%deg)
		

!#######################################################################
!	Diz qual eh a comp gigante
!#######################################################################

		write(*,*) ''
		write(*,*) 'a componente gigante original eh: ', comp_gigante		
		
		write(*,*) ''
				
		write(*,*) 'a componente gigante nao imunizada eh: ', comp_gigante_din

!#######################################################################
!	Diz qual a proporcao da comp gigante antiga estah na nova.
!#######################################################################
		
		write(*,*) 'a proporcao da componente gigante original que nao esta imunizada eh: ', 1.0_dp * comp_gigante_din/comp_gigante
		
!#######################################################################
!	Desaloca as listas ja usadas dentro do escopo da subrotina
!#######################################################################	
	
		deallocate(vertices_a_testar)
		deallocate(rede_original)
		deallocate(lista_tam_clusters_din)

!#######################################################################
!	Aloca lista de distribuicao de graus
!#######################################################################

		
		if(allocated(pok)) deallocate(pok)
			allocate(pok(this%nodes))
			pok = 0.0_dp
			
		if(allocated(grau2)) deallocate(grau2)
			allocate(grau2(this%nodes))
			grau2 = 0
			
!#######################################################################
!	Calculamos distribuicao de grau dentro da comp gigante
!#######################################################################
		
		do i1 = 1, this%nodes
			
			if(lista_de_clusters_din(i1) /= i_comp_gigante_din) cycle
			
			do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
			
				if(sigma(this%listAdj(i2)) /= 2) grau2(i1) = grau2(i1) + 1
			
			enddo
		enddo	

		grauMax2 = maxval(grau2)	
!#######################################################################
!	Calculamos distribuicao p(k).
!#######################################################################
		
		do i1 = 1, size(pok)
			if(grau2(i1) > 0) pok(grau2(i1)) = pok(grau2(i1)) + 1.0_dp
		enddo
		
		nos_grau1 = int(pok(1))
		nos_grau2 = int(pok(2))
		
		sumgrau2 = sum(grau2)
		grauMedio2 = 1.0_dp * sumgrau2/comp_gigante_din
		
        this%degMean = grauMedio2
                
		write(*,*) 'O grau medio efetivo da rede e: ', grauMedio2
		
		write(sorteado_char, '(I0)') sorteado
		
		open(113, file='distribuicao_grau_efetiva_din_'//trim(adjustl(sorteado_char))//'.dat', status = 'unknown')
		
		sorteado = sorteado + 1
!#######################################################################
!	Escreve p(k) em um arquivo.
!#######################################################################
		
		do i1 = 1, size(pok)
			if(pok(i1) > 0.0_dp)then
				pok(i1) = pok(i1)/comp_gigante_din
				write(113,*) i1, pok(i1)
			endif
		enddo

		grauMin2 = minval(grau2)
		if(grauMin2 == 0)then
			do i1=1,size(pok)
				if(pok(i1) > 0.0_dp)then
					grauMin2 = i1
					exit
				endif
			enddo
		endif
		close(113)


!#######################################################################

		grauMax = 1
		do i1 = 1, this%nodes
			if( lista_de_clusters_din(i1) == i_comp_gigante_din)then
				if( grau2(i1) > grauMax)then
					grauMax = grau2(i1)
				endif 
			endif
		enddo

!#######################################################################
!	Diz qual o novo grau maximo
!#######################################################################
		
		write(*,*) 'O grau maximo e: ', grauMax


		
	end subroutine sub_classifica_clusters_nao_imune

!###########################################################################################################################################
	subroutine calcula_k_nn_din(this,criaArquivo, label, arquivo)
	class(grafo) :: this
	integer, intent(in) :: label
	character(len=*) :: arquivo
	real(dp), allocatable :: k_hist(:), k_nn(:)
	real(dp) :: ki_aux, ki_aux2, ki_aux3
	real(dp) :: hist_tam
	logical :: criaArquivo
	integer :: i, j, k, ki_min, ki_max
	
	ki_min = minval(this%deg)
	ki_max = maxval(this%deg)
	
	allocate(k_hist(ki_min:ki_max))
	
	k_hist = 0
	
	do i = 1, this%nodes
		if(lista_de_clusters_din(i) /= i_comp_gigante_din) cycle
		k_hist(grau2(i)) = k_hist(grau2(i)) + 1
	enddo
	
!	if(criaArquivo) then
!		call abreArquivo(label, "k_hist.dat")
!		do i = ki_min, ki_max
!			write(label,*) i, k_hist(i)
!		enddo
!		close(20)
!	endif
	
	hist_tam = sum(k_hist)

	allocate(k_nn(grauMin2:grauMax2))
	k_nn = 0.0_dp
	
	do i = 1, this%nodes
		
		if(lista_de_clusters(i) /= i_comp_gigante) cycle								!Calcula correlacao para cada no i
		
		ki_aux = 0.0_dp
		
		do j = this%aux(i), this%aux(i) + this%deg(i) - 1
			if(sigma(this%listAdj(j)) == 2) cycle
			ki_aux2 = 1d0 * grau2(this%listAdj(j))
			ki_aux = ki_aux + ki_aux2
		enddo
		
		if(grau2(i) /= 0) then
			ki_aux3 = grau2(i)
			k_nn(grau2(i)) = k_nn(grau2(i)) + ki_aux / ki_aux3
		endif
		
		k_hist(grau2(i)) = k_hist(grau2(i)) + 1

	enddo
	
	
		open(label, file=arquivo, status='unknown')
	
		do i = grauMin2, grauMax2
			if(k_hist(i) /= 0) then
				k_nn(i) = 1d0 * k_nn(i) / k_hist(i)
				write(label,*) i, k_nn(i)
			endif
		enddo
		close(label)
				
	deallocate(k_hist)
	
	end subroutine
	

!#######################################################################

	subroutine clustering_din(this, criaArquivo, label, arquivo)

		class(grafo) :: this
		integer, intent(in) :: label
		character(len=*):: arquivo
		real(dp) :: aux, cluster, cluster_global
		integer :: ki_aux
		integer, allocatable :: k_hist(:), C_k(:)
		logical :: criaArquivo
		integer :: i, j, k, l, ki_min, ki_max
		
	cluster = 0_dp
	
	allocate(k_hist(minval(grau2):maxval(grau2)))
	
	k_hist = 0
	
	do i = 1, this%nodes
		if(lista_de_clusters_din(i) /= i_comp_gigante_din) cycle
		k_hist(grau2(i)) = k_hist(grau2(i)) + 1
	enddo
	
	allocate(C_k(minval(grau2):maxval(grau2)))
	
	C_k = 0_dp
	
	do i = 1, this%nodes
		
		if(lista_de_clusters_din(1) /= i_comp_gigante_din) cycle

		if(grau2(i) < 2) cycle
			
				ki_aux = grau2(i)
				aux = 0_dp
				do j = this%aux(i), this%aux(i) + this%deg(i) - 1
					
					if(sigma(this%listAdj(j)) == 2) cycle
					
					do k = this%aux(this%listAdj(j)), this%aux(this%listAdj(j)) + this%deg(this%listAdj(j)) - 1
						
						if(sigma(this%listAdj(k)) == 2) cycle						
						
						if( this%listAdj(k) == i) cycle
						
						do l = j + 1, this%aux(i) + this%deg(i) - 1
							
							if(sigma(this%listAdj(l)) == 2) cycle
							
							if(this%listAdj(l) == this%listAdj(k)) then
								aux = aux + 1
							endif
						enddo
						
					enddo
				enddo

				C_k(grau2(i)) = C_k(grau2(i)) + 2d0 * aux /(ki_aux * (ki_aux - 1))
				cluster = cluster + 2d0 * aux /(ki_aux * (ki_aux - 1))
	enddo
		
		open(label, file=arquivo, status='unknown')
	
		do i = minval(grau2), maxval(grau2)

			if(lista_de_clusters_din(i) /= i_comp_gigante_din) cycle

			if(k_hist(i) /= 0) then
				if(C_k(i) /= 0 )then
					C_k(i) = 1.0_dp * C_k(i) / k_hist(i)
					write(label, *) i, C_k(i)
				endif
			endif
		enddo
	
		close(label)
	
	!	write(*,*) "O valor maximo de list_adj e: ", maxval(list_adj)
	!	write(*,*) "O valor maximo de pos_list e: ", maxval(pos_list)
	!	write(*,*) "O tamanho de list_adj e: ", size(list_adj)


	cluster_global = 1d0 * cluster / comp_gigante_din
			
	write(*,*) "O coeficiente de clusteamento global efetivo eh: ", cluster_global	
	
	deallocate(k_hist)
	
	end subroutine
	
	
end module
