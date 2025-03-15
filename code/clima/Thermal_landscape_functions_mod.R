tdt.curve <- function(ta,time){ 
  
  data <- data.frame(ta,time)
  data <- data[order(data$ta,data$time),]
  #ordena dades per tractaments (de menys a m?s), i dins de tractament, per temps de mort
  
  # Step 1: Calculate CTmax and z from TDT curve
  ta <- as.numeric(levels(as.factor(data$ta)))			
  model <- lm(log10(data$time) ~ data$ta)
  ctmax <- -coef(model)[1]/coef(model)[2]
  z <- -1/coef(model)[2]
  list(model = model, ctmax = ctmax, z = z)}




# Function to estimate thermal tolerance landscape from static assays
# General procedure
# Step 1: Calculate CTmax and z from TDT curve
# Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
# Step 3: Interpolating survival probabilities to make them comparable across treatments 
# Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
# Step 5: Build expected survival curve with mean x and y pooling all data
# Step 6: Expand expected curve to multiple Ta (with 0.1ºC difference for predictive purposes)

#ta: vector amb el tractament t?rmic aplicat a cada individu
#time: vector amb el temps de mort de cada individu

	tolerance.landscape <- function(ta,time){ 
	
			data <- data.frame(ta,time)
			data <- data[order(data$ta,data$time),]
			#ordena dades per tractaments (de menys a m?s), i dins de tractament, per temps de mort

			# Step 1: Calculate CTmax and z from TDT curve
			ta <- as.numeric(levels(as.factor(data$ta)))			
			model <- lm(log10(data$time) ~ data$ta); summary(model)
			ctmax <- -coef(model)[1]/coef(model)[2]
			z <- -1/coef(model)[2]

			# Step 2: Calculate average log10 time and Ta (mean x and y for interpolation purposes)
			time.mn <- mean(log10(data$time))
			ta.mn <- mean(data$ta)
	
			# Step 3: Interpolating survival probabilities to make them comparable across treatments 
			time.interpol <- matrix(,1001,length(ta))
			#matriu de 0s, de 1001 files i una columna per cada tractament
			for(i in 1:length(ta)){	
			time <- c(0,sort(data$time[data$ta==ta[i]])); p <- seq(0,100,length.out = length(time))
			#estimem una probabilitat de superviv?ncia descendent a cada nou individu que va morint
			time.interpol[,i] <- approx(p,time,n = 1001)$y}
			#amb els punts que tinguem time vs p, interpolem linearment el temps de mort que tocaria a cada descens de superviv?ncia de 0.001
			#es tracta de S(time) de cada Ta

			# Step 4: Overlap all survival curves into a single one by shifting each curve to mean x and y employing z
			# Step 5: Build expected survival curve with median survival time for each survival probability
			shift <- (10^((ta - ta.mn)/z)) #vector
			time.interpol.shift <- t(t(time.interpol)*shift)[-1,] #per multiplicar, per cada columna, 1a fila amb 1a fila, 2a fila amb 2a, etc. --> el shift que correspon a cada tractament
			surv.pred <- 10^apply(log10(time.interpol.shift),1,median) 	
	
			# Step 6: Expand predicted survival curves to measured Ta (matrix m arranged from lower to higher ta)
			# Step 7: Obtain predicted values comparable to each empirical measurement
			m <- surv.pred*matrix ((10^((ta.mn - rep(ta, each = 1000))/z)), nrow = 1000) #no ?s multiplicaci? de matrius. fila 1 x fila 1 de cada columna. Predicci? per a cada tractament
			out <-0
			for(i in 1:length(ta)){
				time <- c(0,data$time[data$ta==ta[i]]); p <- seq(0,100,length.out = length(time))
				out <- c(out,approx(seq(0,100,length.out = 1000),m[,i],xout=p[-1])$y)}
				data$time.pred <- out[-1]
				colnames(m) <- paste("time.at",ta,sep=".")
				m <- cbind(surv.prob=seq(1,0.001,-0.001),m)

			par(mfrow=c(1,2),mar=c(4.5,4,1,1),cex.axis=1.1)
			plot(-10,-10,las=1,xlab="Time (min)",ylab="Survival (%)",col="white",xaxs="i",yaxs="i",xlim=c(0,max(data$time)*1.05),ylim=c(0,105))
				for(i in 1:length(ta)){
				time <- c(0,sort(data$time[data$ta==ta[i]])); p <- seq(100,0,length.out = length(time))
				points(time,p,pch=21,bg="black",cex=0.5)
				time <- c(0,sort(data$time.pred[data$ta==ta[i]]))
				points(m[,i+1],100*m[,1],type="l",lty=2)}
				segments(max(data$time)*0.7,90,max(data$time)*0.8,90,lty=2)
				text(max(data$time)*0.82,90,"fitted",adj=c(0,0.5))
			plot(log10(data$time.pred),log10(data$time),pch=21,bg="black",cex=0.5,lwd=0.7,las=1,xlab="Fitted Log10 time",ylab="Measured Log10 time")
			abline(0,1,lty=2)
			rsq <- round(summary(lm(log10(data$time) ~ log10(data$time.pred)))$r.square,3)
			text(min(log10(data$time.pred)),max(log10(data$time)),substitute("r"^2*" = "*rsq),adj=c(0,1))
			list(ctmax = as.numeric(ctmax), z = as.numeric(z), ta.mn = ta.mn,  S = data.frame(surv=seq(0.999,0,-0.001),time=surv.pred),
			time.obs.pred=cbind(data$time,data$time.pred), rsq = rsq)}					

 

# Function to estimate survival probability from tolerance landscapes and environmental temperature data


	dynamic.landscape.mod <- function(ta, tolerance.landscape, sens, sp, day, plot){
			surv <- tolerance.landscape$S[,2]
			ta.mn <- tolerance.landscape$ta.mn
			z <- tolerance.landscape$z
			shift <- 10^((ta.mn - ta)/z)	
			time.rel <- 0
			alive <- 100
			for(i in 1:length(ta)){			
				if(alive[length(na.omit(alive))] > 0) {							
					alive[i+1] <- try(approx(c(0,shift[i]*surv),seq(100,0,length.out = length(c(0,surv))),xout = time.rel[i] + 1)$y,silent=TRUE)		
					time.rel[i+1] <- try(approx(seq(100,0,length.out = length(c(0,surv))),c(0,shift[i + 1]*surv),xout = alive[i + 1])$y,silent=TRUE)}
				else{
					alive <- 0}}				
			out <- tibble(ta = ta,time = 1:length(ta), alive = alive[-1]) 
			#tindria més sentit quedar-se amb la darrera dada d'alive, no la primera, no?
			#equivalent als supervivents abans que passi el minut d'exposició o després.
			
			if (plot) {
			  pdf(file = paste0("figures/daily_mortality/daily_curve_", sens, sp, day, ".pdf"), onefile = T)
			  par(mar=c(4,4,1,1),mfrow=c(1,2))
			  plot(1:length(ta),ta,type="l",xlim=c(0,length(ta)),ylim=c(min(ta),max(ta)),col="black",lwd=1.5,las=1,
				    xlab = "Time (min)", ylab = "Temperature (ºC)")
			  plot(out$time,out$alive,type="l",xlim=c(0,length(ta)),ylim=c(0,100),col="black",lwd=1.5,las=1,
				    xlab = "Time (min)", ylab = "Survival (%)")
			  dev.off()
			}
			return(out)
			}

					

# Function to predict time based on dynamic curve to compare against observed data


	predict.dynamic <- function(time.obs,dynamic.landscape){
			time <- c(1,sort(time.obs)); p <- seq(100,0,length.out = length(time))
			time.pred <- approx(dynamic.landscape$alive,dynamic.landscape$time,xout=p)$y
			list(time.obs.pred = cbind(time[-1],time.pred[-1]))}

