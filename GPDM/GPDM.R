## 1. Carrega os pacotes necessários
library(R2jags)
require(rstan) # apenas para o cálculo de Rhat 

## 2. Define algumas funções auxiliares
# rho.func : retorna uma matriz com os hiperparâmetros rho
rho.func = function(Q,Hj){
	Cj = 2^rowSums(Q)
	center.point = 1
	J = nrow(Q)
	Cmax = max(Cj)
	Hmax = max(Hj)
	rho = array(NA, dim = c(J,Cmax,Hmax))
	for(jj in 1:J){
		cj = Cj[jj]
		hj = Hj[jj]
		cat1 = cat2 = cat3 = vector(length = hj)
		cat1[1] = cj
		cat2[1] = (cj+1)/2
		cat3[1] = 1
		rsm = (cj+1) - c(cat1[1],cat2[1], cat3[1])
		cat1[2:hj] = rep(rsm[1],hj-1)/(hj-1)
		cat2[2:hj] = rep(rsm[2],hj-1)/(hj-1)
		cat3[2:hj] = rep(rsm[3],hj-1)/(hj-1)		
		cinter = setdiff(1:cj,c(1,cj))  # obtain the midway categories vector
		rho[jj,1,1:h1] = cat1
		for(cc in cinter){
			rho[jj,cc,1:h1] = cat2
		}
		rho[jj,cj,1:h1] = cat3
	}
	return(rho)
}
# getmode: retorna a moda de um vetor
getmode = function(v) {
	uniqv = unique(v)
	m =uniqv[which.max(tabulate(match(v, uniqv)))]
	return(m)
}
# patterns.generator: retorna todos os vetores de atributos latentes
patterns.generator=function(K){
	K = as.integer(K)
	C=2^K
	M=matrix(0,nrow=C, ncol=K)
	initial.pos = 1
	for(k in 1:K){
		index.comb=combn(K,k)
		for(ii in 1:ncol(index.comb)){
			pos = initial.pos + ii
			M[pos,index.comb[,ii]]=1
		}
		initial.pos = ii+initial.pos
	}
	colnames(M)=paste0("A",1:K)
	rownames(M)=1:C
	return(M)
}
# reduced.alpha.equivalency: retorna, para cada item, o grupo latente l  
# correspondente a cada classe c
reduced.alpha.equivalency=function(Q)
{
	J=nrow(Q)
	K=ncol(Q)
	C=2^K
	Kj=rowSums(Q)
	Cj=2^Kj
	ID=matrix(0,nrow=J,ncol=C)
	full.alpha.vec=patterns.generator(K)
	for(j in 1:J){
		index=which(Q[j,]==1)
		latent.groups=patterns.generator(Kj[j])
		t.latent.groups=apply(latent.groups,1,paste0,collapse="")
		reduced.alpha = full.alpha.vec[,index,drop=F]
		t.reduced.alpha=apply(reduced.alpha,1,paste0,collapse="")
		ID[j,]=match(t.reduced.alpha,t.latent.groups)
	}
	return(ID)
}
#class.to.vec: recebe um vetor com o índice das classes latentes dos indivíduos 
# e retorna a matriz de vetores de atributos.
class.to.vec=function(pop,K){
	N=length(pop)
	M=matrix(nrow=N,ncol=K)
	alpha_patterns= patterns.generator(K)
	M=t(sapply(pop,function(x) alpha_patterns[x,]))
	colnames(M)=paste0("A",1:K)
	return(M)
}
# prob.matrix.convert.poly: transforma a matriz de probabilidade das classes 
# latentes em uma matriz com os parametros theta.
prob.matrix.convert.poly=function(prob,Q){
	J=nrow(Q)
	K=ncol(Q)
	maxcat = dim(prob)[3]
	C=2^K
	ID=reduced.alpha.equivalency(Q)
	Kj=rowSums(Q)
	Cj=2^Kj
	theta=matrix(NA,nrow = J,ncol=C)
	for(j in 1:J){
		for(c in 1:Cj[j]){
			for(r in 1:maxcat){
				index=match(c, ID[j,])
				theta[j, c, r]=prob[j, index, r]
			}
		}
	}
	return(theta)
}
## 3. Lendo os dados 
# Salve a matriz de resposta no arquivo "Data.txt", a matriz Q no arquivo
# "Q-matrix.txt" e o modelo JAGS no arquivo "GPDM.jags".
# Coloque esses arquivos em uma mesma pasta e copie e 
# cole o caminho dessa pasta na variável load.folder.
#Obs. As respostas devem ser formatadas como números inteiros 
# iniciando do 1.
load.folder = "C:/"  

Q=read.table(file.path(load.folder,"Q-matrix.txt"), sep=",",header = T)
Y=read.table(file.path(load.folder,"Data.txt"), sep=",",header = T)
Q = as.matrix(Q)
Y = as.matrix(Y)
N = nrow(Y)
J = nrow(Q)
K = ncol(Q)
C=2^K
Model=file.path(load.folder,"GPDM.jags")
# Definindo algumas quantidades auxiliares
alpha_patterns= patterns.generator(K)
Lambda = rep(1, C)
Alpha.Equivalency=reduced.alpha.equivalency(Q)
Cj=2^(rowSums(Q))
# Hj é um vetor de dimensao J onde Hj[j] e o numero de niveis 
# de respostas do item j. 
# Forneca Hj ou deixe que seja obtido a partir dos dados 
# atraves do comando da linha abaixo
Hj  = apply(Y,2,  max)
maxcat = max(Hj)
## 4. Ajustando o modelo	
# MCMC setup
n.Iter =  27000
n.Burn =  7000
n.Chains = 2
Thin = 1
jags.inits = NULL
jags.data = list("N"=N, "J"=J,"C"= C,"maxcat"=maxcat,"Hj" =Hj ,"Y"=Y,"Lambda"=Lambda,"Cj"=Cj,"Alpha.Equivalency"=Alpha.Equivalency,"phi"= phi)
jags.parameters = c("prob","ppi","alpha")
GPDMfit = jags( data = jags.data, 
inits = jags.inits,
parameters.to.save = jags.parameters,
model.file =Model,
n.chains = n.Chains, 
n.iter = n.Iter, 
n.burnin=n.Burn,
n.thin = Thin,
DIC = T) 
## 5. Obtendo as estimativas dos parametros
# obtem a classe latente estimada dos individuos
alpha.samples=GPDMfit$BUGSoutput$sims.list$alpha
est.alpha.class=apply(alpha.samples,2,getmode)  
# obtem o vetor de atributos estimados dos individuos 	
est.alpha.vec =class.to.vec(est.alpha.class,K)
# obtem o vetor de proporcoes
est.ppi = GPDMfit$BUGSoutput$mean$ppi
# obtem as estimativas das probabilidades das respostas
# de cada classe latente em cada item
est.prob = GPDMfit$BUGSoutput$mean$prob
patt = patterns.generator(K)
patt = apply(patt,1,paste0,collapse="")
columns.names = paste0("P(",patt,")")
dimnames(est.prob)[1] = paste("item",1:J)
dimnames(est.prob)[2] = columns.names
dimnames(est.prob)[3] = paste("response",1:maxcat)
#obtem os parametros dos itens theta
est.theta = prob.matrix.convert.poly(est.prob,Q)
# DIC
GDINA.DIC =GPDMfit$BUGSoutput$DIC
#obtem o valor de Rhat
sampled.values = GPDMfit$BUGSoutput$sims.array
est.Rhat = apply(sampled.values, 3, rstan::Rhat)
est.Rhat = data.frame(parameter =names(est.Rhat), 
Rhat = unname(est.Rhat))
## 6. salvando as estimativas em arquivos.
# cria uma pasta
save.folder = file.path(load.folder,"Estimates")
if(!dir.exists(save.folder))dir.create(save.folder)
#salva as estimativas
write.csv(est.alpha.class, file=file.path(save.folder,"est.alpha.class.csv"))
write.csv(est.prob, file=file.path(save.folder,"est.prob.csv"))
write.csv(est.theta, file=file.path(save.folder,"est.theta.csv"))
write.csv(est.ppi, file=file.path(save.folder,"est.ppi.csv"))
write(GDINA.DIC, file=file.path(save.folder,"DIC.txt"))
write.csv(est.Rhat, file=file.path(save.folder,"est.Rhat.csv"))
