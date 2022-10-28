module mod_tools
    use geraRede
    implicit none
	
	!###################################################################
	!	Parametros da rede	
	!###################################################################
	integer :: nos_grau1, nos_grau2	
		
		
	!#####################################
	!	Estatistica de tubos
	!#####################################
	integer, private, parameter :: dp = kind(0.0d0)
	real(dp), allocatable :: tamTubos(:)

	!########################################################################
	!	Lista de nos com o tamanho do tubo onde eles fazem parte
	!########################################################################
	
	integer, allocatable :: no_vs_tamanho_tubo(:), hubs(:)
	real(dp) ::	tam_medio_tubo
	real(dp) :: prob_acu_tubos
	!#####################################
	!	Estatistica de clusters
	!#####################################
	integer, allocatable :: lista_de_clusters(:)
	integer :: i_comp_gigante, i_comp_gigante_din, comp_gigante_orig, comp_gigante, comp_gigante_din
	real(dp) :: prop_cg_din

	!#####################################
	!	Estatistica de distancias
	!#####################################

	real(dp) :: dist_media, dist_media_k, k_corte
	integer :: maior_dist
    integer, allocatable :: lista_distancias(:)

    contains


!###################################################################################
!		Estatistica de tubos
!###################################################################################

        subroutine statistic_Tubes(this)


            !#######################################################################
            !	Objeto do tipo grafo
            !#######################################################################	
            class(grafo) :: this

            !#######################################################################
            !	Parametro de precisao
            !#######################################################################	

 
            !#######################################################################
            !	Variaveis mudas
            !#######################################################################	
            integer :: i1, i2, i3, i4

            !#######################################################################
            !	Variaveis internas
            !#######################################################################	
            integer :: degDif2


            !#######################################################################
            !	Lista degDif2
            !#######################################################################	
            integer, allocatable :: l_degDif2(:), mask(:)

            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: pDT(:)


            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: auxPDT(:)

            !#######################################################################
            !	Lista que guarda o numero de pontas de tubo que sai de um no
            !	degDif2. Sempre que fizermos uma busca, vamos checar
            !	pDT ateh nPDT
            !
            !#######################################################################	

            integer, allocatable :: nPDT(:)
			

            !#######################################################################
            !	Lista que guarda o numero atual de tubos registrados
            !	que saem de um no degDif2
            !#######################################################################	
            integer :: ultimoTubo
            integer :: sumDegDif2
            integer :: buffer, vizDaVez

            !#######################################################################
            !	Histograma tamanho de tubos
            !#######################################################################	
			integer :: normaTamTubos
            integer :: tamTuboProv
            logical :: novo

            !######################################################
            !	Contamos o numero de nos com deg /= 2
            !	Aqui estou contando mesmo os nos de grau 1
            !######################################################	

            degDif2 = 0
            sumDegDif2 = 0
            do i1 = 1, this%nodes
	            if(this%deg(i1) /= 2)then
					degDif2 = degDif2 + 1
					sumDegDif2 = sumDegDif2 + this%deg(i1)
				endif
            enddo
			
		
            !######################################################
            !	Aloca lista de degDif2
            !######################################################
            if(allocated(l_degDif2)) deallocate(l_degDif2)
	            allocate(l_degDif2(degDif2))


            !######################################################
            !	Aloca lista de pontas de tubo encontradas
            !	nas vizinhancas dos nos Dif2
            !	seu acesso sera semelhante aa this%listAdj,
            !	feito por meio de this%aux.
            !######################################################
            !	Mas ela precisa da lista rede%aux para acessar,
            !	no entanto, nao podemos usar junto a rede%deg,
            !	portanto, temos que inventar outra, com a mesma
            !	funcao da rede%deg, mas que comece zerada.
            !	Ela vai guardar o numero de vizinhos de grau
            !	2 registrados, afim de que numa proxima busca,
            !	nao comecemos a procura por tubos atraves
            !	de um stub repetido e que levara a um tubo ja
            !	contado.
            !######################################################



            if(allocated(pDT)) deallocate(pDT)
	            allocate(pDT(sumDegDif2))


            if(allocated(nPDT)) deallocate(nPDT)
	            allocate(nPDT(degDif2))
	            nPDT = 0


            if(allocated(auxPDT)) deallocate(auxPDT)
	            allocate(auxPDT(degDif2))

			if(allocated(mask)) deallocate(mask)
				allocate(mask(this%nodes))

            !######################################################
            !	Zero degDif2 pra recontar e colocar
            !	nos degDif2 na lista
            !######################################################


            degDif2 = 0
            do i1 = 1, this%nodes
	            if(this%deg(i1) /= 2) then
		            degDif2 = degDif2 + 1
		            mask(i1) = degDif2
		            l_degDif2(degDif2) = i1		
	            endif
            enddo	




            auxPDT(1) = 1
            do i1 = 2, degDif2
	            auxPDT(i1) = auxPDT(i1 - 1) + this%deg(l_degDif2(i1 - 1))
            enddo


            !######################################################
            !	Lista cujos indices sao os tamanhos dos tubos
            !	E cujas entradas sao o numero de vezes
            !	Que cada tamanho aparece
            !######################################################

            if(allocated(tamTubos)) deallocate(tamTubos)
	            allocate(tamTubos(this%nodes - degDif2))
	            tamTubos = 0.0_dp





            !#######################################################################
            !	Aqui comeca o algoritmo
            !#######################################################################	

            !##############################################################
            !	Olhamos todos os nos que nao formam tubo
            !##############################################################

            do i1 = 1, degDif2

	            !#############################################################
	            !	Procuramos pontas novas de tubos dentre seus
	            !	vizinhos.
	            !#############################################################
	            
            olhaViz:	do i2 = this%aux(l_degDif2(i1)), this%aux(l_degDif2(i1)) + this%deg(l_degDif2(i1)) - 1

		            
		            
			            !#################################################
			            !	Se ele tem grau dois, pode ser um tubo.
			            !#################################################		
					            
			            if(this%deg(this%listAdj(i2)) /= 2) cycle olhaViz
				            
		            
		            
		            
				            !#################################################
				            ! Damos um voto de confianca
				            !#################################################

				            novo = .True.
				            
				            
				            !#################################################
				            ! Mas checamos se ele eh final de tubo.
				            !#################################################
				            ! Elaborar a ideia da auxPDT foi dificil,
				            ! apesar de a ideia vir pronta com a ipos do
				            ! Wesley. Entender que ela nao era a mesma coisa
				            ! que a ipos ou a minha this%aux demorou um pouco.
				            !#################################################
				            
				            do i3 = 1, nPDT(i1)
					            if(this%listAdj(i2) == pDT(auxPDT(i1) + i3 - 1)) novo = .False.
				            enddo
				            
				            !############################################################
				            ! Se ele continua novo, somamos um
				            ! no tamanho do tubo.
				            !############################################################

				            
				            
				            if(novo)then 
					            tamTuboProv = 1
				            else
					            cycle olhaViz
				            endif
				            
				            !############################################################
				            ! Agora temos que olhar o vizinho do 
				            ! tubo tam 1 para saber se ele pode ser 2,
				            ! ou 3, ou ...
				            !############################################################
				            ! No proximo loop, temos que evitar o no
				            ! de origem e refluxo.
				            ! Se o no da vez eh o de origem, faz nada e segue
				            ! o rumo.
				            ! 
				            ! Se eh um no ja visitado, guardado no
				            ! buffer, nao faz nada e segue o rumo olhando
				            ! outro. Se eh um no degDif2 diferente do de
				            ! origem, guarda o tamTuboProv
				            ! na lista tamTubo,
				            ! tamTubo(tamTuboProv) = tamTubo(tamTuboProv) + 1,
				            ! pDT(auxPDT(l_degDif2(i1)) + nPDT(l_degDif2(i1))) = degDif2
				            ! faz nPDT(l_degDif2(i1)) = nPDT(l_degDif2(i1)) + 1
				            ! e exit contaMissangas.
				            !############################################################
				            
				            buffer = l_degDif2(i1)
				            vizDaVez = this%listAdj(i2)

            contaMissangas2:	do
            contaMissangas:		do i3 = this%aux(vizDaVez), this%aux(vizDaVez) + this%deg(vizDaVez) - 1
						            !###########################################################
						            ! Vamos olhar os vizinhos de this%listAdj(i2)
						            !###########################################################
						            if(this%listAdj(i3) /= buffer)then
							            if(this%deg(this%listAdj(i3)) == 2)then
								            tamTuboProv = tamTuboProv + 1
								            buffer = vizDaVez
								            vizDaVez = this%listAdj(i3)
								            exit contaMissangas
							            else
								            tamTubos(tamTuboProv) = tamTubos(tamTuboProv) + 1
								            pDT(auxPDT(mask(this%listAdj(i3))) + nPDT(mask(this%listAdj(i3)))) = vizDaVez
								            nPDT(mask(this%listAdj(i3))) = nPDT(mask(this%listAdj(i3))) + 1
								            exit contaMissangas2
							            endif
						            endif
					            enddo contaMissangas
				            enddo contaMissangas2
	            enddo olhaViz	
            enddo
            
            normaTamTubos = sum(tamTubos)
            
            do i3 = 1, size(tamTubos)
				tamTubos(i3) = 1.0_dp * tamTubos(i3)/normaTamTubos
           enddo
        end subroutine


!#######################################################################
        subroutine statistic_Tubes_c_detalhe(this, label, arquivo)


            !#######################################################################
            !	Objeto do tipo grafo
            !#######################################################################	
            class(grafo) :: this

			integer :: label
			character(len=300) :: arquivo
            !#######################################################################
            !	Parametro de precisao
            !#######################################################################	

 
            !#######################################################################
            !	Variaveis mudas
            !#######################################################################	
            integer :: i1, i2, i3, i4, i5, pos_no

            !#######################################################################
            !	Variaveis internas
            !#######################################################################	
            integer :: degDif2


            !#######################################################################
            !	Lista degDif2
            !#######################################################################	
            integer, allocatable :: l_degDif2(:), mask(:)

            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: pDT(:)


            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: auxPDT(:)

            !#######################################################################
            !	Lista que guarda o numero de pontas de tubo que sai de um no
            !	degDif2. Sempre que fizermos uma busca, vamos checar
            !	pDT ateh nPDT
            !
            !#######################################################################	

            integer, allocatable :: nPDT(:)
			
			!########################################################################
			!	Lista provisoria que vai guardar os nos de grau dois que foram 
			!	encontrados no tubo atual.
			!########################################################################
			
			integer, allocatable :: nos_tubo_atual(:)
			
			
            !#######################################################################
            !	Lista que guarda o numero atual de tubos registrados
            !	que saem de um no degDif2
            !#######################################################################	
            integer :: ultimoTubo
            integer :: sumDegDif2
            integer :: buffer, vizDaVez

            !#######################################################################
            !	Histograma tamanho de tubos
            !#######################################################################	
			integer :: normaTamTubos
            integer :: tamTuboProv
            logical :: novo

            !######################################################
            !	Contamos o numero de nos com deg /= 2
            !	Aqui estou contando mesmo os nos de grau 1
            !######################################################	

            degDif2 = 0
            sumDegDif2 = 0
            do i1 = 1, this%nodes

				if(lista_de_clusters(i1) /= i_comp_gigante) cycle
	    
	            if(this%deg(i1) == 2) cycle

				degDif2 = degDif2 + 1

				sumDegDif2 = sumDegDif2 + this%deg(i1)
			
            enddo
			
			
			!###########################################################
			!	Essa eh a lista que vai guardar os nos com grau 2
			!	e seu tamanho de tubo.
			!	Uma forma de fazer isso eh, quando for imunizar os nos
			!	de grau dois em funcao do tamanho, primeiro detecto
			!	seu grau. Se ele tem grau 2, olho nessa lista por
			!	seu tamanho no local que corresponde ao seu indice.
			!###########################################################
			
			
			if(allocated(no_vs_tamanho_tubo)) deallocate(no_vs_tamanho_tubo)
				allocate(no_vs_tamanho_tubo(this%nodes))
		
			if(allocated(nos_tubo_atual)) deallocate(nos_tubo_atual)
				allocate(nos_tubo_atual(this%nodes))
		
            !######################################################
            !	Aloca lista de degDif2
            !######################################################
            if(allocated(l_degDif2)) deallocate(l_degDif2)
	            allocate(l_degDif2(degDif2))


            !######################################################
            !	Aloca lista de pontas de tubo encontradas
            !	nas vizinhancas dos nos Dif2
            !	seu acesso sera semelhante aa this%listAdj,
            !	feito por meio de this%aux.
            !######################################################
            !	Mas ela precisa da lista rede%aux para acessar,
            !	no entanto, nao podemos usar junto a rede%deg,
            !	portanto, temos que inventar outra, com a mesma
            !	funcao da rede%deg, mas que comece zerada.
            !	Ela vai guardar o numero de vizinhos de grau
            !	2 registrados, afim de que numa proxima busca,
            !	nao comecemos a procura por tubos atraves
            !	de um stub repetido e que levara a um tubo ja
            !	contado.
            !######################################################



            if(allocated(pDT)) deallocate(pDT)
	            allocate(pDT(sumDegDif2))


            if(allocated(nPDT)) deallocate(nPDT)
	            allocate(nPDT(degDif2))
	            nPDT = 0


            if(allocated(auxPDT)) deallocate(auxPDT)
	            allocate(auxPDT(degDif2))

			if(allocated(mask)) deallocate(mask)
				allocate(mask(this%nodes))
				mask = 0
            !######################################################
            !	Zero degDif2 pra recontar e colocar
            !	nos degDif2 na lista
            !######################################################


            degDif2 = 0
            do i1 = 1, this%nodes
		
				if(lista_de_clusters(i1) /= i_comp_gigante) cycle
	    
	            if(this%deg(i1) == 2) cycle
		
				degDif2 = degDif2 + 1
		
		        mask(i1) = degDif2
		
		        l_degDif2(degDif2) = i1
            enddo	




            auxPDT(1) = 1
            do i1 = 2, degDif2
	            auxPDT(i1) = auxPDT(i1 - 1) + this%deg(l_degDif2(i1 - 1))
            enddo


            !######################################################
            !	Lista cujos indices sao os tamanhos dos tubos
            !	E cujas entradas sao o numero de vezes
            !	Que cada tamanho aparece
            !######################################################

            if(allocated(tamTubos)) deallocate(tamTubos)
	            allocate(tamTubos(this%nodes - degDif2))
	            tamTubos = 0.0_dp





            !#######################################################################
            !	Aqui comeca o algoritmo
            !#######################################################################	

            !##############################################################
            !	Olhamos todos os nos que nao formam tubo
            !##############################################################

            do i1 = 1, degDif2

	            !#############################################################
	            !	Procuramos pontas novas de tubos dentre seus
	            !	vizinhos.
	            !#############################################################
	            
            olhaViz:	do i2 = this%aux(l_degDif2(i1)), this%aux(l_degDif2(i1)) + this%deg(l_degDif2(i1)) - 1

		            
		            
			            !#################################################
			            !	Se ele tem grau dois, pode ser um tubo.
			            !#################################################		
					            
			            if(this%deg(this%listAdj(i2)) /= 2) cycle olhaViz
				            
		            
		            
		            
				            !#################################################
				            ! Damos um voto de confianca
				            !#################################################

				            novo = .True.
				            
				            
				            !#################################################
				            ! Mas checamos se ele eh final de tubo.
				            !#################################################
				            ! Elaborar a ideia da auxPDT foi dificil,
				            ! apesar de a ideia vir pronta com a ipos do
				            ! Wesley. Entender que ela nao era a mesma coisa
				            ! que a ipos ou a minha this%aux demorou um pouco.
				            !#################################################
				            
				            do i3 = 1, nPDT(i1)
					            if(this%listAdj(i2) == pDT(auxPDT(i1) + i3 - 1)) novo = .False.
				            enddo
				            
				            !############################################################
				            ! Se ele continua novo, somamos um
				            ! no tamanho do tubo.
				            !############################################################

				            
				            
				            if(novo)then 
					            tamTuboProv = 1
				            else
					            cycle olhaViz
				            endif
				            
				            !############################################################
				            ! Agora temos que olhar o vizinho do 
				            ! tubo tam 1 para saber se ele pode ser 2,
				            ! ou 3, ou ...
				            !############################################################
				            ! No proximo loop, temos que evitar o no
				            ! de origem e refluxo.
				            ! Se o no da vez eh o de origem, faz nada e segue
				            ! o rumo.
				            ! 
				            ! Se eh um no ja visitado, guardado no
				            ! buffer, nao faz nada e segue o rumo olhando
				            ! outro. Se eh um no degDif2 diferente do de
				            ! origem, guarda o tamTuboProv
				            ! na lista tamTubo,
				            ! tamTubo(tamTuboProv) = tamTubo(tamTuboProv) + 1,
				            ! pDT(auxPDT(l_degDif2(i1)) + nPDT(l_degDif2(i1))) = degDif2
				            ! faz nPDT(l_degDif2(i1)) = nPDT(l_degDif2(i1)) + 1
				            ! e exit contaMissangas.
				            !############################################################
				            
				            buffer = l_degDif2(i1)
				            vizDaVez = this%listAdj(i2)
							
							!########################################################
								pos_no = 1
								nos_tubo_atual(pos_no) = this%listAdj(i2)
							!########################################################
									
            contaMissangas2:	do
            contaMissangas:		do i3 = this%aux(vizDaVez), this%aux(vizDaVez) + this%deg(vizDaVez) - 1
						            !###########################################################
						            ! Vamos olhar os vizinhos de this%listAdj(i2)
						            !###########################################################
						            if(this%listAdj(i3) /= buffer)then
							            if(this%deg(this%listAdj(i3)) == 2)then
								            tamTuboProv = tamTuboProv + 1
								            buffer = vizDaVez
								            vizDaVez = this%listAdj(i3)
								            
								            !#########################################
												pos_no = pos_no + 1
												nos_tubo_atual(pos_no) = this%listAdj(i3)
								            !#########################################
							
								            exit contaMissangas
							            else
								            tamTubos(tamTuboProv) = tamTubos(tamTuboProv) + 1
								            pDT(auxPDT(mask(this%listAdj(i3))) + nPDT(mask(this%listAdj(i3)))) = vizDaVez
								            nPDT(mask(this%listAdj(i3))) = nPDT(mask(this%listAdj(i3))) + 1
								            exit contaMissangas2
							            endif
						            endif
					            enddo contaMissangas
				            enddo contaMissangas2
				            
							!#########################################
								do i5 = 1, pos_no
									no_vs_tamanho_tubo(nos_tubo_atual(i5)) = tamTuboProv
								enddo
							!#########################################
				            
	            enddo olhaViz	
            enddo
            
            normaTamTubos = sum(tamTubos)
            
            do i3 = 1, size(tamTubos)
				tamTubos(i3) = 1.0_dp * tamTubos(i3)/normaTamTubos
			enddo
           
           
			arquivo = trim(adjustl(arquivo))
			
           	open(label, file=arquivo,status='unknown')
			
			do i3 = 1, size(tamTubos)
				if(tamTubos(i3) > 0.0_dp)then
					write(label,*) i3, tamTubos(i3)
				endif	
			enddo
			
			close(label)
           
           
            !###########################################################
            !	Calcula tubo medio
            !###########################################################
           
			tam_medio_tubo = 0.0_dp
           
			do i3 = 1, size(tamTubos)
				if(tamTubos(i3) > 0.0_dp)then
					tam_medio_tubo = tam_medio_tubo + 1.0_dp * i3 * tamTubos(i3)
				endif
			enddo
           
			write(*,*) " "
			write(*,*) "o tamanho medio dos tubos e: ", tam_medio_tubo
			write(*,*) " "
			
			
			prob_acu_tubos = 0.0_dp	
			
           do i3 = 1, size(tamTubos)
			if( 1.0_dp * i3 < tam_medio_tubo)then
				prob_acu_tubos = prob_acu_tubos + tamTubos(i3)
			else
				exit
			endif
           enddo
           
           write(*,*) " "
           write(*,*) "A probabilidade de encontrar um tubo menor que o tam medio ", tam_medio_tubo, " e ", prob_acu_tubos
           write(*,*) " "
           
           deallocate(nos_tubo_atual)
           
        end subroutine

!#######################################################################
        subroutine statistic_Tubes_n_imun_c_detalhe(this)


            !#######################################################################
            !	Objeto do tipo grafo
            !#######################################################################	
            class(grafo) :: this

            !#######################################################################
            !	Parametro de precisao
            !#######################################################################	

 
            !#######################################################################
            !	Variaveis mudas
            !#######################################################################	
            integer :: i1, i2, i3, i4, i5, pos_no

            !#######################################################################
            !	Variaveis internas
            !#######################################################################	
            integer :: degDif2


            !#######################################################################
            !	Lista degDif2
            !#######################################################################	
            integer, allocatable :: l_degDif2(:), mask(:)

            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: pDT(:)


            !#######################################################################
            !	Lista pra guardar os nos de grau 2 nas vizinhancas de um no Dif2
            !#######################################################################	

            integer, allocatable :: auxPDT(:)

            !#######################################################################
            !	Lista que guarda o numero de pontas de tubo que sai de um no
            !	degDif2. Sempre que fizermos uma busca, vamos checar
            !	pDT ateh nPDT
            !
            !#######################################################################	

            integer, allocatable :: nPDT(:)
			
			!########################################################################
			!	Lista provisoria que vai guardar os nos de grau dois que foram 
			!	encontrados no tubo atual.
			!########################################################################
			
			integer, allocatable :: nos_tubo_atual(:)
			
			
            !#######################################################################
            !	Lista que guarda o numero atual de tubos registrados
            !	que saem de um no degDif2
            !#######################################################################	
            integer :: ultimoTubo
            integer :: sumDegDif2
            integer :: buffer, vizDaVez

            !#######################################################################
            !	Histograma tamanho de tubos
            !#######################################################################	
			
			integer :: normaTamTubos
            integer :: tamTuboProv
            logical :: novo

            !######################################################
            !	Contamos o numero de nos com deg /= 2
            !	Aqui estou contando mesmo os nos de grau 1
            !######################################################	

            degDif2 = 0
       
            sumDegDif2 = 0
       
            do i1 = 1, this%nodes
		
				if(lista_de_clusters_din(i1) /= i_comp_gigante_din) cycle
	    
	            if(this%deg(i1) == 2) cycle

				degDif2 = degDif2 + 1
	
				sumDegDif2 = sumDegDif2 + this%deg(i1)
			
            enddo
			
			
			!###########################################################
			!	Essa eh a lista que vai guardar os nos com grau 2
			!	e seu tamanho de tubo.
			!	Uma forma de fazer isso eh, quando for imunizar os nos
			!	de grau dois em funcao do tamanho, primeiro detecto
			!	seu grau. Se ele tem grau 2, olho nessa lista por
			!	seu tamanho no local que corresponde ao seu indice.
			!###########################################################
			
			
			if(allocated(no_vs_tamanho_tubo)) deallocate(no_vs_tamanho_tubo)
				allocate(no_vs_tamanho_tubo(this%nodes))
		
			if(allocated(nos_tubo_atual)) deallocate(nos_tubo_atual)
				allocate(nos_tubo_atual(this%nodes - degDif2))
		
            !######################################################
            !	Aloca lista de degDif2
            !######################################################
            if(allocated(l_degDif2)) deallocate(l_degDif2)
	            allocate(l_degDif2(degDif2))


            !######################################################
            !	Aloca lista de pontas de tubo encontradas
            !	nas vizinhancas dos nos Dif2
            !	seu acesso sera semelhante aa this%listAdj,
            !	feito por meio de this%aux.
            !######################################################
            !	Mas ela precisa da lista rede%aux para acessar,
            !	no entanto, nao podemos usar junto a rede%deg,
            !	portanto, temos que inventar outra, com a mesma
            !	funcao da rede%deg, mas que comece zerada.
            !	Ela vai guardar o numero de vizinhos de grau
            !	2 registrados, afim de que numa proxima busca,
            !	nao comecemos a procura por tubos atraves
            !	de um stub repetido e que levara a um tubo ja
            !	contado.
            !######################################################



            if(allocated(pDT)) deallocate(pDT)
	            allocate(pDT(sumDegDif2))


            if(allocated(nPDT)) deallocate(nPDT)
	            allocate(nPDT(degDif2))
	            nPDT = 0


            if(allocated(auxPDT)) deallocate(auxPDT)
	            allocate(auxPDT(degDif2))

			if(allocated(mask)) deallocate(mask)
				allocate(mask(this%nodes))

            !######################################################
            !	Zero degDif2 pra recontar e colocar
            !	nos degDif2 na lista
            !######################################################


            degDif2 = 0

            do i1 = 1, this%nodes

				if(lista_de_clusters_din(i1) /= i_comp_gigante_din) cycle
	    
	            if(this%deg(i1) == 2) cycle

				degDif2 = degDif2 + 1

				mask(i1) = degDif2

				l_degDif2(degDif2) = i1		


            enddo	




            auxPDT(1) = 1
            do i1 = 2, degDif2
	            auxPDT(i1) = auxPDT(i1 - 1) + this%deg(l_degDif2(i1 - 1))
            enddo


            !######################################################
            !	Lista cujos indices sao os tamanhos dos tubos
            !	E cujas entradas sao o numero de vezes
            !	Que cada tamanho aparece
            !######################################################

            if(allocated(tamTubos)) deallocate(tamTubos)
	            allocate(tamTubos(this%nodes - degDif2))
	            tamTubos = 0.0_dp





            !#######################################################################
            !	Aqui comeca o algoritmo
            !#######################################################################	

            !##############################################################
            !	Olhamos todos os nos que nao formam tubo
            !##############################################################

            do i1 = 1, degDif2

	            !#############################################################
	            !	Procuramos pontas novas de tubos dentre seus
	            !	vizinhos.
	            !#############################################################
	            
            olhaViz:	do i2 = this%aux(l_degDif2(i1)), this%aux(l_degDif2(i1)) + this%deg(l_degDif2(i1)) - 1
						
						
		            
		            
			            !#################################################
			            !	Se ele tem grau dois, pode ser um tubo.
			            !#################################################		
					            
			            if(this%deg(this%listAdj(i2)) /= 2) cycle olhaViz
				            
		            		            
				            !#################################################
				            ! Damos um voto de confianca
				            !#################################################

				            novo = .True.
				            
				            
				            !#################################################
				            ! Mas checamos se ele eh final de tubo.
				            !#################################################
				            ! Elaborar a ideia da auxPDT foi dificil,
				            ! apesar de a ideia vir pronta com a ipos do
				            ! Wesley. Entender que ela nao era a mesma coisa
				            ! que a ipos ou a minha this%aux demorou um pouco.
				            !#################################################
				            
				            do i3 = 1, nPDT(i1)
					            if(this%listAdj(i2) == pDT(auxPDT(i1) + i3 - 1)) novo = .False.
				            enddo
				            
				            !############################################################
				            ! Se ele continua novo, somamos um
				            ! no tamanho do tubo.
				            !############################################################

				            
				            
				            if(novo)then 
					            tamTuboProv = 1
				            else
					            cycle olhaViz
				            endif
				            
				            !############################################################
				            ! Agora temos que olhar o vizinho do 
				            ! tubo tam 1 para saber se ele pode ser 2,
				            ! ou 3, ou ...
				            !############################################################
				            ! No proximo loop, temos que evitar o no
				            ! de origem e refluxo.
				            ! Se o no da vez eh o de origem, faz nada e segue
				            ! o rumo.
				            ! 
				            ! Se eh um no ja visitado, guardado no
				            ! buffer, nao faz nada e segue o rumo olhando
				            ! outro. Se eh um no degDif2 diferente do de
				            ! origem, guarda o tamTuboProv
				            ! na lista tamTubo,
				            ! tamTubo(tamTuboProv) = tamTubo(tamTuboProv) + 1,
				            ! pDT(auxPDT(l_degDif2(i1)) + nPDT(l_degDif2(i1))) = degDif2
				            ! faz nPDT(l_degDif2(i1)) = nPDT(l_degDif2(i1)) + 1
				            ! e exit contaMissangas.
				            !############################################################
				            
				            buffer = l_degDif2(i1)
				            vizDaVez = this%listAdj(i2)
							
							!########################################################
								pos_no = 1
								nos_tubo_atual(pos_no) = this%listAdj(i2)
							!########################################################
									
            contaMissangas2:	do
            contaMissangas:		do i3 = this%aux(vizDaVez), this%aux(vizDaVez) + this%deg(vizDaVez) - 1
									
									if(sigma(this%listAdj(i3)) == 2) exit contaMissangas2
									
						            !###########################################################
						            ! Vamos olhar os vizinhos de this%listAdj(i2)
						            !###########################################################
						            if(this%listAdj(i3) /= buffer)then
							            if(this%deg(this%listAdj(i3)) == 2)then
								            tamTuboProv = tamTuboProv + 1
								            buffer = vizDaVez
								            vizDaVez = this%listAdj(i3)
								            
								            !#########################################
												pos_no = pos_no + 1
												nos_tubo_atual(pos_no) = this%listAdj(i3)
								            !#########################################
							
								            exit contaMissangas
							            else
								            tamTubos(tamTuboProv) = tamTubos(tamTuboProv) + 1
								            pDT(auxPDT(mask(this%listAdj(i3))) + nPDT(mask(this%listAdj(i3)))) = vizDaVez
								            nPDT(mask(this%listAdj(i3))) = nPDT(mask(this%listAdj(i3))) + 1
								            exit contaMissangas2
							            endif
						            endif
					            enddo contaMissangas
				            enddo contaMissangas2
				            
						!#########################################
							do i5 = 1, pos_no
								no_vs_tamanho_tubo(nos_tubo_atual(i5)) = tamTuboProv
							enddo
						!#########################################
				            
	            enddo olhaViz	
            enddo
            
            normaTamTubos = sum(tamTubos)
            
            do i3 = 1, size(tamTubos)
				tamTubos(i3) = 1.0_dp * tamTubos(i3)/normaTamTubos
			enddo
           
			open(111, file='lista_tam_tubos_n_imun.dat',status='unknown')
			
			do i3 = 1, size(tamTubos)
				if(tamTubos(i3) > 0.0_dp)then
					write(111,*) i3, tamTubos(i3)
				endif
			enddo
			
			close(111)
			
            !###########################################################
            !	Calcula tubo medio
            !###########################################################
           
			tam_medio_tubo = 0.0_dp
           
			do i3 = 1, size(tamTubos)
				if(tamTubos(i3) > 0.0_dp)then
					tam_medio_tubo = tam_medio_tubo + 1.0_dp * i3 * tamTubos(i3)
				endif
			enddo
           
			write(*,*) " "
			write(*,*) "o tamanho medio dos tubos e: ", tam_medio_tubo
			write(*,*) " "
			
			
			prob_acu_tubos = 0.0_dp	
			
           do i3 = 1, size(tamTubos)
			if( 1.0_dp * i3 < tam_medio_tubo)then
				prob_acu_tubos = prob_acu_tubos + tamTubos(i3)
			else
				exit
			endif
           enddo
           
           write(*,*) " "
           write(*,*) "A probabilidade de encontrar um tubo menor que o tam medio ", tam_medio_tubo, " e ", prob_acu_tubos
           write(*,*) " "
           
           deallocate(nos_tubo_atual)
           
        end subroutine

!#######################################################################


!#################################################################################################################################################	
!		Calcula e classifica tamanho das componentes
!#################################################################################################################################################
	
	
	subroutine sub_classifica_clusters(this, criaArquivo, label, arquivo)
	class(grafo), intent(inout) :: this
	integer, intent(in) :: label
	character(len = *) :: arquivo
	logical :: criaArquivo
	real(dp) :: grauMedio2
	integer ::  sumgrau2
	character(len=100) :: sorteado_char
	!############################################################
	!	Variaveis mudas
	!############################################################
		integer :: i1, i2, i3, i4
	!############################################################
	!	Listas auxiliares internas
	!############################################################
		
		integer, allocatable :: lista_tam_clusters(:), vertices_a_testar(:), rede_original(:)

	!############################################################
	!	Variaveis internas
	!############################################################
		integer :: n_zeros, indice_cluster, novidade, agregados, pos_novosVertices, pos_ultimosVertices




!#################################################################################################################################################	
!		Alocamento de listas auxiliares internas
!#################################################################################################################################################


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


	if(allocated(lista_de_clusters)) deallocate(lista_de_clusters)
		allocate(lista_de_clusters(this%nodes))
		lista_de_clusters = 0
						!Aqui vão entrar os índices dos clusters. Acho que nao precisa, mas settei como zero, por seguranca.
						!Se nao tiver dado errado, se eu der uma busca por zeros, nao vai aparecer nada.


	if(allocated(lista_tam_clusters)) deallocate(lista_tam_clusters)
		allocate(lista_tam_clusters(this%nodes))

		lista_tam_clusters = 0



				
		

	!###################################################################################################
	!				Inicializacoes importantes
	!###################################################################################################

			vertices_a_testar(1) = 1

			indice_cluster = 1			!O primeiro cluster tem indice 1.
				
			n_zeros = 0

			pos_novosVertices = 1
			pos_ultimosVertices = 1


			lista_de_clusters(vertices_a_testar(1)) = indice_cluster
				
			rede_original(vertices_a_testar(1)) = 0		
					
			lista_tam_clusters(indice_cluster) = lista_tam_clusters(indice_cluster) + 1
			n_zeros = n_zeros + 1
	

		do while(n_zeros < this%nodes)
		
			novidade = 1				! Necessario inicializar esse cara como nao nulo
		

			do while(novidade /= 0)
						
															 				
					agregados = pos_ultimosVertices
				
					novidade = 0
				
					do i1 = pos_novosVertices, pos_ultimosVertices
						
							do i2 = this%aux(vertices_a_testar(i1)), this%aux(vertices_a_testar(i1)) + this%deg(vertices_a_testar(i1)) - 1
														
								if(rede_original(this%listAdj(i2)) /= 0)then
							
									agregados = agregados + 1
									rede_original(this%listAdj(i2)) = 0
									n_zeros = n_zeros + 1
									vertices_a_testar(agregados) = this%listAdj(i2)
							
									lista_de_clusters(this%listAdj(i2)) = indice_cluster
						
									lista_tam_clusters(indice_cluster) = lista_tam_clusters(indice_cluster) + 1					
											
									novidade = novidade + 1	
						
								endif
							enddo
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
			
					indice_cluster = indice_cluster + 1					!Quer dizer que ha um novo cluster
			
					vertices_a_testar(pos_novosVertices) = i2				! Como nao zero mais, o novo a testar
														! 
					lista_de_clusters(i2) = indice_cluster
				
					rede_original(i2) = 0							! Antes eu zerava esse cara
														!
					
					lista_tam_clusters(indice_cluster) = lista_tam_clusters(indice_cluster) + 1
			
					n_zeros = n_zeros + 1

					exit  loop_p
			
				endif
			enddo loop_p
		enddo
	
	if(criaArquivo) then
		call abreArquivo(label, arquivo)

		do i1 =1, this%nodes
			if(lista_tam_clusters(i1) /= 0) then
				write(label,*) i1, lista_tam_clusters(i1)
			endif
		enddo

			close(label)

			call abreArquivo(label + 1, arquivo)
			write(label + 1, *) "Vertice ", "Cluster"
	
		do i1 = 1, size(lista_de_clusters)
			write(label + 1,*) i1, lista_de_clusters(i1)
		enddo
	
		close(label + 1)
	endif
	
		comp_gigante = maxval(lista_tam_clusters)
		

		do i1 = 1, this%nodes
			if(lista_tam_clusters(i1) == comp_gigante)then

				i_comp_gigante = i1
				
				exit
				
			endif
		enddo
		
		i2 = 0
		do i1 = 1, this%nodes
			if(lista_tam_clusters(i1) == comp_gigante)then
				i2 = i2 + 1	
			endif
			if(i2 > 1) stop "Ha duas componentes gigantes para esta semente. Interrompendo processo."
		enddo
				
		
		write(*,*) 'O tamanho da componente gigante e: ', comp_gigante
		
		
		
		deallocate(vertices_a_testar)
		deallocate(rede_original)
	
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
			if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
					grau2(i1) = grau2(i1) + 1
			enddo
		enddo	



!#######################################################################
!	Calculamos distribuicao p(k).
!#######################################################################
		
		do i1 = 1, size(pok)
			if(grau2(i1) > 0)then
				pok(grau2(i1)) = pok(grau2(i1)) + 1.0_dp
			endif
		enddo
		
		nos_grau1 = int(pok(1))
		nos_grau2 = int(pok(2))
				
		sumgrau2 = sum(grau2)
		grauMedio2 = 1.0_dp * sumgrau2/comp_gigante
		
               this%degMean = grauMedio2
                
		write(*,*) 'O grau medio efetivo da rede e: ', grauMedio2
		write(*,*) ' '
		write(*,*) 'O grau maximo da componente gigante e: ', maxval(grau2)
		
		write(sorteado_char, '(I0)') sorteado
		sorteado_char = trim(adjustl(sorteado_char))
		
		open(113, file='distribuicao_grau_efetiva_'//trim(adjustl(sorteado_char))//'.dat', status = 'unknown')
		
		sorteado = sorteado + 1
!#######################################################################
!	Escreve p(k) em um arquivo.
!#######################################################################
		
		do i1 = 1, size(pok)
			if(pok(i1) > 0.0_dp)then
				pok(i1) = pok(i1)/comp_gigante
				write(113,*) i1, pok(i1)
			endif
		enddo

		close(113)
	
	
	end subroutine sub_classifica_clusters
	


!#################################################################################################################################################

		
	
!#################################################################################################################################################
	
!#################################################################################################################################################	
!		Calcula distancias medias e diametro da rede
!#################################################################################################################################################
	
		subroutine calcula_distancias(this)
		implicit none
		integer, parameter :: dp=kind(0.0d0)
		integer, allocatable :: vertices_a_testar(:)
		integer, allocatable :: lista_distancias(:)
		integer(kind=1), allocatable :: rede_original(:)
		integer :: n_novos, indice_cluster, testados, pos_novosVertices, pos_ultimosVertices, dist, &
		maior_dist_prov, cont_dist, N_dist, N_3_core, N_k_core
		
		

		!###################################################################################################
		!	Obejto da glasse grafo
		!###################################################################################################		
		class(grafo) :: this
		
		!###################################################################################################
		!	variaveis auxiliares
		!###################################################################################################
		integer :: i1, i2, i3

		!###################################################################################################
		!	Alocamento de listas internas
		!	Fazer umas gambiarras aqui... Vou ter que fazer uma lista de nos da componente gigante
		!	apenas...
		!	Olhar a lista de clusters e procurar por nos que tenham o indice da componente gigante.
		!	O tamanho das listas que vem a seguir terao tamanho = comp_gigante, obvio.
		!	Vou alocar as listas com este tamanho
		!###################################################################################################
		
		if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina sub_classifica_clusters, do modulo mod_tool_redes. Finalizando programa."

		if(allocated(lista_distancias)) deallocate(lista_distancias)
			allocate(lista_distancias(this%nodes))
										

	
		if(allocated(rede_original))	deallocate(rede_original)
			allocate(rede_original(this%nodes))
		
										! Se alguem eh elegivel, ou seja, se ainda nao tiver sido visitado, recebe 1. Do contrario, 0.
	
		if(allocated(vertices_a_testar)) deallocate(vertices_a_testar)	
			allocate(vertices_a_testar(this%nodes))
			
	
		!###################################################################################################
		!	Inicializo variaveis importantes
		!###################################################################################################
!		N_dist = comp_gigante * (comp_gigante - 1)
		maior_dist = 0 
		dist_media = 0.0_dp
		dist_media_k = 0.0_dp		
		
	
		N_3_core = 0
		N_k_core = 0
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < 3) cycle
			N_3_core = N_3_core + 1
		enddo
	
		write(*,*) "O tamanho do 3-core e ", N_3_core
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < k_corte) cycle
			N_k_core = N_k_core + 1
		enddo

		write(*,*) "O tamanho do k-core e ", N_k_core

		
		!###################################################################################################
		!	Comeco o algoritmo de busca
		!###################################################################################################
		do i1 = 1, this%nodes
		
			!#########################################################################
			!		Olhar a lista de cluster... se o noh i1 tem indice i_comp_cigante,
			!		segue, senao, cycle. Solucao mais que simples.
			!#########################################################################
			if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			
			!if(this%deg(i1) < 3) cycle
			
			rede_original = 1					! Estou dizendo que todo mundo eh elegivel pra ser visitado	
			rede_original(i1) = 0					! menos o cara por onde eu comeco
			
			pos_ultimosVertices = 1
			pos_novosVertices = 1
			
			vertices_a_testar = 0
			vertices_a_testar(1) = i1
			
			lista_distancias = 0
			
			dist = 0
			n_novos = 1						! Parece que settei 1 aqui soh pra entrar no loop
			
			do while(n_novos > 0)
				
				n_novos = 0
				
				testados = pos_ultimosVertices
				
				dist = dist + 1
				
				do i2 = pos_novosVertices, pos_ultimosVertices
					
					do i3 = this%aux(vertices_a_testar(i2)), this%aux(vertices_a_testar(i2)) + this%deg(vertices_a_testar(i2)) - 1
					
						if(this%deg(this%listAdj(i3)) == 1)then
							rede_original(this%listAdj(i3)) = 0
							cycle
						endif
						
						if(rede_original(this%listAdj(i3)) /=0) then
							
							rede_original(this%listAdj(i3)) = 0
							
							testados = testados + 1				! Aqui dizemos que testamos a distancia de mais
																! um.
							vertices_a_testar(testados) = this%listAdj(i3)
							
							n_novos = n_novos + 1				! Aqui dizemos que um novo noh com 
																! possiveis novos vizinhos foi encontrado.
																															
							lista_distancias(this%listAdj(i3)) = dist
						endif
					enddo
				enddo
				pos_novosVertices = pos_ultimosVertices + 1
				pos_ultimosVertices = testados
			enddo
		
			maior_dist_prov = maxval(lista_distancias)
			
			maior_dist = max(maior_dist_prov, maior_dist)

!			dist_media = dist_media + 1.0_dp * sum(1.0_dp * lista_distancias)/(comp_gigante - 1)

			
			if(this%deg(i1) >= k_corte)then
				do i2 = 1, this%nodes
					if(this%deg(i2) >= k_corte)then
						dist_media_k = dist_media_k + 1.0_dp * lista_distancias(i2) / (N_k_core - 1)
					endif
				enddo
			else
				do i2 = 1, this%nodes
					if(this%deg(i2) > 2)then
						dist_media = dist_media + 1.0_dp * lista_distancias(i2) / (N_3_core - 1)
					endif
				enddo			
			endif
			
		enddo
		
!		dist_media = 1.0_dp * dist_media / comp_gigante
		
		dist_media = 1.0_dp * dist_media / N_3_core
		
		dist_media_k = 1.0_dp * dist_media_k / N_k_core
				
		write(*,*) 'A distancia media da componente gigante e: ', dist_media
		write(*,*) ' '

		write(*,*) 'A distancia media entre hubs e: ', dist_media_k
		write(*,*) ' '
		write(*,*) 'Onde os hubs sao tomados quando o grau for maior que ', k_corte
		write(*,*) ' '


		write(*,*) 'O diametro da componente gigante e: ', maior_dist
		write(*,*) ' '	
		
		deallocate(rede_original)
		deallocate(vertices_a_testar)
			
	end subroutine calcula_distancias


!=======================================================================
		subroutine calcula_distancias_ao_hub(this)
		implicit none
		integer, parameter :: dp=kind(0.0d0)
		integer, allocatable :: vertices_a_testar(:)
		integer(kind=1), allocatable :: rede_original(:)
		integer :: n_novos, indice_cluster, testados, pos_novosVertices, pos_ultimosVertices, dist, &
		maior_dist_prov, cont_dist, N_dist, N_3_core, N_k_core
		
		

		!###################################################################################################
		!	Obejto da glasse grafo
		!###################################################################################################		
		class(grafo) :: this
		
		!###################################################################################################
		!	variaveis auxiliares
		!###################################################################################################
		integer :: i1, i2, i3

		!###################################################################################################
		!	Alocamento de listas internas
		!	Fazer umas gambiarras aqui... Vou ter que fazer uma lista de nos da componente gigante
		!	apenas...
		!	Olhar a lista de clusters e procurar por nos que tenham o indice da componente gigante.
		!	O tamanho das listas que vem a seguir terao tamanho = comp_gigante, obvio.
		!	Vou alocar as listas com este tamanho
		!###################################################################################################
		
		if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina sub_classifica_clusters, do modulo mod_tool_redes. Finalizando programa."

		if(allocated(lista_distancias)) deallocate(lista_distancias)
			allocate(lista_distancias(this%nodes))
										
		if(allocated(rede_original))	deallocate(rede_original)
			allocate(rede_original(this%nodes))
										! Se alguem eh elegivel, ou seja, se ainda nao tiver sido visitado, recebe 1. Do contrario, 0.
		if(allocated(vertices_a_testar)) deallocate(vertices_a_testar)	
			allocate(vertices_a_testar(this%nodes))
	
		!###################################################################################################
		!	Inicializo variaveis importantes
		!###################################################################################################
!		N_dist = comp_gigante * (comp_gigante - 1)
		maior_dist = 0 
		dist_media = 0.0_dp
		dist_media_k = 0.0_dp		
		
	
		N_3_core = 0
		N_k_core = 0
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < 3) cycle
			N_3_core = N_3_core + 1
		enddo
	
		write(*,*) "O tamanho do 3-core e ", N_3_core
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < k_corte) cycle
			N_k_core = N_k_core + 1
		enddo

		write(*,*) "O tamanho do k-core e ", N_k_core

		
		!###################################################################################################
		!	Comeco o algoritmo de busca
		!###################################################################################################
		!---------------------------------------------------------------
		!  Esse eh o hub.
		!---------------------------------------------------------------
		i1 = 1
		
		!if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			
		rede_original = 1					! Estou dizendo que todo mundo eh elegivel pra ser visitado	
		rede_original(i1) = 0					! menos o cara por onde eu comeco
		
		pos_ultimosVertices = 1
		pos_novosVertices = 1
		
		vertices_a_testar = 0
		vertices_a_testar(1) = i1
		
		lista_distancias = 0

		dist = 0
		n_novos = 1						! Parece que settei 1 aqui soh pra entrar no loop

		do while(n_novos > 0)
		
			n_novos = 0
		
			testados = pos_ultimosVertices

			dist = dist + 1

			do i2 = pos_novosVertices, pos_ultimosVertices
				do i3 = this%aux(vertices_a_testar(i2)), this%aux(vertices_a_testar(i2)) + this%deg(vertices_a_testar(i2)) - 1

					if(rede_original(this%listAdj(i3)) == 0) cycle

					if(this%deg(this%listAdj(i3)) == 1)then
						rede_original(this%listAdj(i3)) = 0
						lista_distancias(this%listAdj(i3)) = dist
						cycle
					endif

					rede_original(this%listAdj(i3)) = 0

					testados = testados + 1				! Aqui dizemos que testamos a distancia de mais
															! um.
					vertices_a_testar(testados) = this%listAdj(i3)

					n_novos = n_novos + 1				! Aqui dizemos que um novo noh com 
														! possiveis novos vizinhos foi encontrado.

					lista_distancias(this%listAdj(i3)) = dist
				enddo
			enddo
			!-----------------------------------------------------------
			! Depois que saimos desse loop, significa que testamos
			! desde o pos_novosVertices ateh o pos_ultimosVertices
			! antigos.
			! Assim, o proximo 'pos_novosVertices' tem que ser a partir
			! do pos_ultimosVertices antigos.
			!-----------------------------------------------------------
			pos_novosVertices = pos_ultimosVertices + 1
			pos_ultimosVertices = testados
		enddo

		if(this%deg(i1) >= k_corte)then
				do i2 = 1, this%nodes
					if(this%deg(i2) >= k_corte)then
						dist_media_k = dist_media_k + 1.0_dp * lista_distancias(i2) / (N_k_core - 1)
					endif
				enddo
			else
				do i2 = 1, this%nodes
					if(this%deg(i2) > 2)then
						dist_media = dist_media + 1.0_dp * lista_distancias(i2) / (N_3_core - 1)
					endif
				enddo			
			endif


		maior_dist = maxval(lista_distancias)

!		dist_media = 1.0_dp * dist_media / comp_gigante
		
		dist_media = 1.0_dp * dist_media / N_3_core
		
		dist_media_k = 1.0_dp * dist_media_k / N_k_core
				
		write(*,*) 'A distancia media da componente gigante e: ', dist_media
		write(*,*) ' '

		write(*,*) 'A distancia media entre hubs e: ', dist_media_k
		write(*,*) ' '
		write(*,*) 'Onde os hubs sao tomados quando o grau for maior que ', k_corte
		write(*,*) ' '


		write(*,*) 'O diametro da componente gigante e: ', maior_dist
		write(*,*) ' '	
		
		deallocate(rede_original)
		deallocate(vertices_a_testar)
			
	end subroutine
!=======================================================================
!#################################################################################################################################################
!	Distancias dentro do nucleo SF
		subroutine calcula_distancias_no_nucleo(this)
		implicit none
		integer, parameter :: dp=kind(0.0d0)
		integer, allocatable :: vertices_a_testar(:)
		integer, allocatable :: lista_distancias(:)
		integer(kind=1), allocatable :: rede_original(:)
		integer :: n_novos, indice_cluster, testados, pos_novosVertices, pos_ultimosVertices, dist, &
		maior_dist_prov, cont_dist, N_dist, N_3_core, N_k_core

		

		!###################################################################################################
		!	Obejto da glasse grafo
		!###################################################################################################		
		class(grafo) :: this
		
		!###################################################################################################
		!	variaveis auxiliares
		!###################################################################################################
		integer :: i1, i2, i3
		
		


		!###################################################################################################
		!	Alocamento de listas internas
		!	Fazer umas gambiarras aqui... Vou ter que fazer uma lista de nos da componente gigante
		!	apenas...
		!	Olhar a lista de clusters e procurar por nos que tenham o indice da componente gigante.
		!	O tamanho das listas que vem a seguir terao tamanho = comp_gigante, obvio.
		!	Vou alocar as listas com este tamanho
		!###################################################################################################
		
		if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina sub_classifica_clusters, do modulo mod_tool_redes. Finalizando programa."

		if(allocated(lista_distancias)) deallocate(lista_distancias)
			allocate(lista_distancias(this%nodes))
										

	
		if(allocated(rede_original))	deallocate(rede_original)
			allocate(rede_original(this%nodes))
		
										! Se alguem eh elegivel, ou seja, se ja nao tiver sido visitado, recebe 1. Do contrario, 0.
	
		if(allocated(vertices_a_testar)) deallocate(vertices_a_testar)	
			allocate(vertices_a_testar(this%nodes))
			
	
		!###################################################################################################
		!	Inicializo variaveis importantes
		!###################################################################################################
!		N_dist = comp_gigante * (comp_gigante - 1)
		maior_dist = 0 
		dist_media = 0.0_dp
		dist_media_k = 0.0_dp		
	
		N_3_core = 0
		N_k_core = 0
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < 3) cycle
			N_3_core = N_3_core + 1
		enddo
		
		
		do i1 =1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante)cycle
			if(this%deg(i1) < k_corte) cycle
			N_k_core = N_k_core + 1
		enddo
			
		!###################################################################################################
		!	Comeco o algoritmo de busca
		!###################################################################################################
		do i1 = 1, this%nodes
		
			!#########################################################################
			!		Olhar a lista de cluster... se o noh i1 tem indice i_comp_cigante,
			!		segue, senao, cycle. Solucao mais que simples.
			!#########################################################################
			
			if(this%deg(i1) < 3) cycle
			
			if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			
			rede_original = 1					! Estou dizendo que todo mundo eh elegivel pra ser visitado	
			rede_original(i1) = 0					! menos o cara por onde eu comeco
			
			pos_ultimosVertices = 1					
			pos_novosVertices = 1
			
			vertices_a_testar = 0
			vertices_a_testar(1) = i1
			
			lista_distancias = 0
			
			dist = 0
			n_novos = 1						! Parece que settei 1 aqui soh pra entrar no loop
			
			do while(n_novos > 0)
				
				n_novos = 0
				
				testados = pos_ultimosVertices
				
				dist = dist + 1
				
				do i2 = pos_novosVertices, pos_ultimosVertices
					
					do i3 = this%aux(vertices_a_testar(i2)), this%aux(vertices_a_testar(i2)) + this%deg(vertices_a_testar(i2)) - 1
						
						
						if(lista_de_clusters(this%listAdj(i3)) /= i_comp_gigante) cycle
						
						if(this%deg(this%listAdj(i3)) == 1)then
							rede_original(this%listAdj(i3)) = 0
							cycle
						endif
						
						if(rede_original(this%listAdj(i3)) /=0) then
							
							rede_original(this%listAdj(i3)) = 0
							
							testados = testados + 1				! Aqui dizemos que testamos a distancia de mais
																! um.
							vertices_a_testar(testados) = this%listAdj(i3)
							
							n_novos = n_novos + 1				! Aqui dizemos que um novo noh com 
																! possiveis novos vizinhos foi encontrado.
																															
							lista_distancias(this%listAdj(i3)) = dist
						endif
					enddo
				enddo
				pos_novosVertices = pos_ultimosVertices + 1
				pos_ultimosVertices = testados
			enddo
		
			maior_dist_prov = maxval(lista_distancias)
			
			maior_dist = max(maior_dist_prov, maior_dist)

			if(this%deg(i1) >= k_corte)then
				do i2 = 1, size(hubs)
						dist_media_k = dist_media_k + 1.0_dp * lista_distancias(hubs(i2)) / (N_k_core - 1)
				enddo
			else
				do i2 = 1, this%nodes
					if(this%deg(i2) > 2)then
						dist_media = dist_media + 1.0_dp * lista_distancias(i2) / (N_3_core - 1)
					endif
				enddo			
			endif
			
		enddo
		
!		dist_media = 1.0_dp * dist_media / comp_gigante
		
		dist_media = 1.0_dp * dist_media / N_3_core
		
		dist_media_k = 1.0_dp * dist_media_k / N_k_core
				
		write(*,*) 'A distancia media dentro do 3-core e: ', dist_media
		write(*,*) ' '

		write(*,*) 'A distancia media entre hubs e: ', dist_media_k
		write(*,*) ' '
		write(*,*) 'Onde hubs são considerados com grau acima de ', k_corte
		write(*,*) ' '

		write(*,*) 'O diametro do nucleo scale free e: ', maior_dist
		write(*,*) ' '	
		
		deallocate(rede_original)
		deallocate(vertices_a_testar)
			
	end subroutine calcula_distancias_no_nucleo
	
!####################################################################################################################################################		
	
	
	
	!#################################################################################################################################################	
!		Calcula distancias medias e diametro da rede
!#################################################################################################################################################
	
		subroutine calcula_distancias_sampling(this, sampling_prop)
		implicit none
		integer, parameter :: dp=kind(0.0d0)

	
		!###################################################################################################
		!	Argumentos
		!###################################################################################################		
		class(grafo), intent(in) :: this
		real(dp), intent(in) :: sampling_prop
		
		!###################################################################################################
		!	variaveis auxiliares
		!###################################################################################################		
		integer, allocatable :: vertices_a_testar(:), lista_distancias(:)
		integer(kind=1), allocatable :: rede_original(:)
		integer :: n_novos, indice_cluster, testados, pos_novosVertices, pos_ultimosVertices, dist, &
		maior_dist_prov, cont_dist, N_dist
		integer :: sampling_size
		integer :: i1, i2, i3, i4
		integer :: label
		!###################################################################################################
		!	Alocamento de listas internas
		!	Fazer umas gambiarras aqui... Vou ter que fazer uma lista de nos da componente gigante
		!	apenas...
		!	Olhar a lista de clusters e procurar por nos que tenham o indice da componente gigante.
		!	O tamanho das listas que vem a seguir terao tamanho = comp_gigante, obvio.
		!	Vou alocar as listas com este tamanho
		!###################################################################################################
		
		if(.not. allocated(lista_de_clusters)) stop "Voce precisa primeiro rodar a subrotina 'sub_classifica_clusters', do modulo mod_tool_redes.f90. Finalizando programa."

		if(allocated(lista_distancias)) deallocate(lista_distancias)
			allocate(lista_distancias(this%nodes))
										

	
		if(allocated(rede_original)) deallocate(rede_original)
			allocate(rede_original(this%nodes))
		
										! Se alguem eh elegivel, ou seja, se ja nao tiver sido visitado, recebe 1. Do contrario, 0.
	
		if(allocated(vertices_a_testar)) deallocate(vertices_a_testar)	
			allocate(vertices_a_testar(comp_gigante))
			
	
		!#########################################################################################################
		!	Inicializo variaveis importantes
		!#########################################################################################################

		!#########################################################################################################
		! Como a rotina que caracteriza a componente gigante jah foi implementada, podemos usar suas variaveis
		!#########################################################################################################
		
		sampling_size = int(sampling_prop * comp_gigante)
		
				
		maior_dist = 0 
		dist_media = 0.0_dp
	
		!###################################################################################################
		!	Comeco o algoritmo de busca
		!###################################################################################################

		!##################################################################################################################
		!	Inicializo a variavel que vai contar o numero de amostras tiradas da rede criada, ateh completar sampling_size
		!##################################################################################################################
		
		i4 = 0
		
		do i1 = 1, this%nodes
		
			if(i4 >= sampling_size) exit
			
			!#########################################################################
			!		Olhar a lista de cluster... se o noh i1 tem indice i_comp_cigante,
			!		segue, senao, cycle. Solucao mais que simples.
			!#########################################################################
			
			if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			
			
			rede_original = 1					! Estou dizendo que todo mundo eh elegivel pra ser visitado	
			rede_original(i1) = 0					! menos o cara por onde eu comeco
			
			pos_ultimosVertices = 1					
			pos_novosVertices = 1
			
			vertices_a_testar = 0
			vertices_a_testar(1) = i1
			
			lista_distancias = 0
			
			dist = 0
			n_novos = 1						! Parece que settei 1 aqui soh pra entrar no loop
			
			do while(n_novos > 0)
				
				n_novos = 0
				
				testados = pos_ultimosVertices
				
				dist = dist + 1
				
				do i2 = pos_novosVertices, pos_ultimosVertices
					
					do i3 = this%aux(vertices_a_testar(i2)), this%aux(vertices_a_testar(i2)) + this%deg(vertices_a_testar(i2)) - 1
						
						if(rede_original(this%listAdj(i3)) /=0) then
							
							rede_original(this%listAdj(i3)) = 0
							
							testados = testados + 1				! Aqui dizemos que testamos a distancia de mais
																! um.
							vertices_a_testar(testados) = this%listAdj(i3)
							
							n_novos = n_novos + 1				! Aqui dizemos que um novo noh com 
																! possiveis novos vizinhos foi encontrado.
																															
							lista_distancias(this%listAdj(i3)) = dist
						endif
					enddo
				enddo
				pos_novosVertices = pos_ultimosVertices + 1
				pos_ultimosVertices = testados
			enddo
		
			maior_dist_prov = maxval(lista_distancias)
			
			maior_dist = max(maior_dist_prov, maior_dist)

			dist_media = dist_media + 1.0_dp * sum(1.0_dp * lista_distancias)/(comp_gigante - 1)
			
			i4 = i4 + 1
			
		enddo
		dist_media = 1.0_dp * dist_media / sampling_size	
		deallocate(rede_original)
		deallocate(vertices_a_testar)
	end subroutine calcula_distancias_sampling
	
	!###################################################
	!		Abre arquivos
	!###################################################
	
	subroutine abreArquivo(label, nomeArquivo)
		character(len = 255) :: workDir
		character(len = *), intent(in) :: nomeArquivo
		character(len = 510) :: caminho
		integer, intent(in) :: label
		
		call getcwd(workDir)
		caminho = trim(adjustl(trim(adjustl(workDir))//"/"//trim(adjustl(nomeArquivo))))
		open(label, file = caminho, status = "unknown")
	end subroutine abreArquivo
	
!#######################################################################
	
	subroutine clustering(this, criaArquivo, label, arquivo)

		class(grafo) :: this
		integer, intent(in) :: label
		character(len=*):: arquivo
		real(dp) :: aux, cluster, cluster_global
		integer :: ki_aux
		integer, allocatable :: k_hist(:)
		real(dp), allocatable ::  C_k(:)
		logical :: criaArquivo
		integer :: i, j, k, l, ki_min, ki_max
		
	cluster = 0_dp
	
	allocate(k_hist(this%degMin:this%degMax))
	
	k_hist = 0
	
	
	allocate(C_k(this%degMin:this%degMax))
	
	C_k = 0_dp
	
	do i = 1, this%nodes
		
		if(lista_de_clusters(1) /= i_comp_gigante) cycle

		if(this%deg(i) < 2) cycle
			
				ki_aux = this%deg(i)
				aux = 0_dp
				do j = this%aux(i), this%aux(i) + this%deg(i) - 1
				
					do k = this%aux(this%listAdj(j)), this%aux(this%listAdj(j)) + this%deg(this%listAdj(j)) - 1
						if( this%listAdj(k) /= i) then
								do l = j + 1, this%aux(i) + this%deg(i) - 1
									if(this%listAdj(l) == this%listAdj(k)) then
										aux = aux + 1
									endif
								enddo
						endif
					enddo
				enddo
				C_k(this%deg(i)) = C_k(this%deg(i)) + 2_dp * aux /(ki_aux * (ki_aux - 1))
				cluster = cluster + 2_dp * aux /(ki_aux * (ki_aux - 1))
				
				k_hist(this%deg(i)) = k_hist(this%deg(i)) + 1
	enddo
		
		if(criaArquivo)then
			open(label, file=arquivo, status='unknown')
		
			do i = this%degMin, this%degMax

				if(k_hist(i) > 0) then
					if(C_k(i) > 0._dp )then
						C_k(i) = 1._dp * C_k(i) / k_hist(i)
						write(label, *) i, C_k(i)
					endif
				endif
			enddo
		
			close(label)
		endif
		
	!	write(*,*) "O valor maximo de list_adj e: ", maxval(list_adj)
	!	write(*,*) "O valor maximo de pos_list e: ", maxval(pos_list)
	!	write(*,*) "O tamanho de list_adj e: ", size(list_adj)

	cluster_global = 1d0 * cluster / comp_gigante
			
	write(*,*) "O coeficiente de clusteamento global eh: ", cluster_global	
	
	deallocate(k_hist)
	
	end subroutine
	
	
!#######################################################################

	subroutine calcula_k_nn(this,criaArquivo, label, arquivo)
	class(grafo) :: this
	integer, intent(in) :: label
	character(len=*) :: arquivo
	integer, allocatable :: k_hist(:)
	real(dp), allocatable :: k_nn(:)
	real(dp) :: ki_aux, ki_aux2, ki_aux3
	real(dp) :: hist_tam
	logical :: criaArquivo
	integer :: i, j, k, ki_min, ki_max
	
	allocate(k_hist(this%degMin:this%degMax))
	
	k_hist = 0
	
	allocate(k_nn(this%degMin:this%degMax))
	k_nn = 0.0_dp
	
	do i = 1, this%nodes
		
		if(lista_de_clusters(i) /= i_comp_gigante) cycle								!Calcula correlacao para cada no i
		
		ki_aux = 0.0_dp
		
		do j = this%aux(i), this%aux(i) + this%deg(i) - 1
			ki_aux2 = 1d0 * this%deg(this%listAdj(j))
			ki_aux = ki_aux + ki_aux2
		enddo
		
		if(this%deg(i) /= 0) then
			ki_aux3 = this%deg(i)
			k_nn(this%deg(i)) = k_nn(this%deg(i)) + ki_aux / ki_aux3
		endif
		
		k_hist(this%deg(i)) = k_hist(this%deg(i)) + 1

	enddo
		
	hist_tam = sum(k_hist)
		
		
	if(criaArquivo)then	
		open(unit=label, file=arquivo, status='unknown')	
		
		do i = this%degMin, this%degMax
		
			if(k_hist(i) /= 0) then
				k_nn(i) = 1d0 * k_nn(i) / k_hist(i)
				write(label,*) i, k_nn(i)
			endif
		
		enddo
		
		close(label)
	endif
	
	deallocate(k_hist)
	
	end subroutine	
	
end module
