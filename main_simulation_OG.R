# Set source file path
source("../src/source_Experiments.R")

# Set hyperparameter
seed <- sample(x=1:5000,size=1000)
rate <- c(0.001, 0.05, 0.1)
N <- 1000

# Simulation 1
t_knlg=c(1:1000)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:100){

      dt0 <- generate_simulation(simulation = "simulation1", n = 1000, p = 100, outlier_rate = outlier_rate)#, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data
      
      seq2.1 <- seq_transform2(dts = dt, sequence = c("O"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate, save_data =T)
      #seq2.2 <- seq_transform2(dts = seq2.1$transformed_data$O, sequence = c("TK"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate, save_data =T)
      #seq2 <- seq_transform2(dts = seq2.2$transformed_data$TK, sequence = c("DK"),d_knlg=-1, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate, save_data =T)
      
      if(0){
        plotting(dt,outl=outlying)
        plotting( seq2.1$transformed_data$O,outl=outlying)
        plotting( seq2.2$transformed_data$TK,outl=outlying)
        plotting( seq2$transformed_data$DK,outl=outlying)
      }

      seq2 <- seq_transform2(dts = dt, sequence = c("T0"),t_knlg=t_knlg, depth_method = "erld", direction = "both", outlier_rate = outlier_rate)
      ST_rst<-evaluation(seq2$outliers$T0,outlying,N=N)
      
      seq2 <- seq_transform2(dts = dt, sequence = c("T0","DK"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
      DST_rst<-evaluation(seq2$outliers$T0,outlying,N=N)
      
      if(0){
        seq2 <- seq_transform2(dts = dt, sequence = c("DK"),t_knlg=t_knlg,d_knlg=-1, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate)
        DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)
        
        seq2.2 <- seq_transform2(dts = seq2.1$transformed_data$O, sequence = c("DK"),t_knlg=t_knlg,d_knlg=-1, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate)
        DST_rst2<-evaluation(seq2.2$outliers$DK,outlying,N=N)
        
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "T2"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
        
        seq2 <- seq_transform2(dts = dt, sequence = c("DK"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
        DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)
      }

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
      #evl_table<-data.frame(rbind(evl_table,c(msplot_rst,TVD_rst,MUOD_rst,ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("ST_pc","ST_pf","DST_pc","DST_pf" ))#"msplot_pc","msplot_pf","TVD_pc","TVD_pf","MUOD_pc","MUOD_pf","ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}

# Simulation 2
t_knlg=c(1:100)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation2", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data
      
      if(0){
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "T2"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
      }
      seq0 <- seq_transform2(dts = dt, sequence = c("T0"),t_knlg=t_knlg, depth_method = "erld", direction = "both", outlier_rate = outlier_rate)
      ST_rst<-evaluation(seq0$outliers$T0,outlying,N=N)

      seq2 <- seq_transform2(dts = dt, sequence = c("DK"),t_knlg=t_knlg, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate)
      DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
      #evl_table<-data.frame(rbind(evl_table,c(msplot_rst,TVD_rst,MUOD_rst,ST_rst,DST_rst)))
    }
    #names(evl_table)<-paste0(c("msplot_pc","msplot_pf","TVD_pc","TVD_pf","MUOD_pc","MUOD_pf","ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}

# Simulation 3
t_knlg=c(20:30)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation3", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data

      if(0){
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "T2"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
      }
      

      seq0 <- seq_transform2(dts = dt, sequence = c("T1", "T3"),t_knlg=t_knlg, depth_method = "erld", direction = "both", outlier_rate = outlier_rate)
      ST_rst<-evaluation(seq0$outliers$T3,outlying,N=N)

      seq2 <- seq_transform2(dts = dt, sequence = c("T1", "TK", "T3", "DK"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
      DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}

# Simulation 4
t_knlg=c(20:30)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation4", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data

      if(0){
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "T2"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
      }
      
      seq0 <- seq_transform2(dts = dt, sequence = c("T1","T3"),t_knlg=t_knlg, depth_method = "erld", direction = "both", outlier_rate = outlier_rate)
      ST_rst<-evaluation(seq0$outliers$T3,outlying,N=N)

      seq2 <- seq_transform2(dts = dt, sequence = c("T1", "TK","T3", "DK"),t_knlg=t_knlg, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate)
      DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}

# Simulation 5
t_knlg=c(20:40)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation5", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data

      if(0){

        
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "D1"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
      }
      
      seq0 <- seq_transform2(dts = dt, sequence = c("T1","T3"),t_knlg=t_knlg, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate,save=F)
      ST_rst<-evaluation(seq0$outliers$DK,outlying,N=N)
      
      #plotting(seq0$transformed_data$DK,outl=outlying)
      #plotting(seq2$transformed_data$T1,outl=outlying)
      #plotting(dt,outl=outlying)
      
      seq2 <- seq_transform2(dts = dt, sequence = c("T1","TK","T3", "DK"),t_knlg=t_knlg, depth_method = "erld", direction = "lower", outlier_rate = outlier_rate,save=F)
      DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}

# Simulation 6
t_knlg=c(1:100)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation6", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])

      outlying <-dt0$true_outliers
      dt <- dt0$data
      
      if(0){
        msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
        msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)
        
        TVD_object <- tvdmss2(dt)
        TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)
        
        MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
        MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)
        
        seq0 <- seq_transform(dts = dt, sequence = c("T1", "D1"),  depth_method = "erld", erld_type = "two_sided")
        ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)
        
        #plotting(seq0$transformed_data$TK,outl=outlying)
        #plotting(seq2$transformed_data$DK,outl=outlying)
        png(file="step1.png",width=600,height = 600)
        plotting(dt,outl=outlying,xl="Index",yl="Value")
        dev.off()
        png("step2.png",width=600,height = 600)
        plotting(seq2$transformed_data$D1,outl=outlying,xl="Index",yl="Value")
        dev.off()
        png("step3.png",width=600,height = 600)
        plotting(seq2$transformed_data$T5,outl=outlying,xl="Frequency",yl="Value")
        dev.off()
        png("step4.png",width=600,height = 600)
        plotting(seq2$transformed_data$TK,outl=outlying,xaxis=c(2:3), xl="Frequency",yl="Value")
        dev.off()
        
      }
      

      seq0 <- seq_transform2(dts = dt, sequence = c("D1","T5"),t_knlg=c(3:4), depth_method = "erld", direction = "both", outlier_rate = outlier_rate,save=F)
      ST_rst<-evaluation(seq0$outliers$T5,outlying,N=N)
      

      seq2 <- seq_transform2(dts = dt, sequence = c("D1","T5","TK"),t_knlg=c(3:4), depth_method = "erld", direction = "upper", outlier_rate = outlier_rate,save=T)
      DST_rst<-evaluation(seq2$outliers$TK,outlying,N=N) #"TK","T1","O","DK"

      evl_table<-data.frame(rbind(evl_table,c(ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}




# Simulation 6 with structural knowledge
t_knlg=c(1:100)
for (i in rate){
    outlier_rate <- i
    evl_table <- c()
    for (j in 1:1000){

      dt0 <- generate_simulation(simulation = "simulation6", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[j])
      dt1 <- generate_simulation(simulation = "auxiliary2", n = 1000, p = 100, outlier_rate = outlier_rate, seed = seed[i])

      outlying <-dt0$true_outliers
      dt <- dt0$data
      dt1 <- dt1$data

      msplot_object <- msplot(dts = dt, return_mvdir=FALSE, plot=F)
      msplot_rst<-evaluation(msplot_object$outliers,outlying,N=N)

      TVD_object <- tvdmss2(dt)
      TVD_rst<-evaluation(TVD_object$outliers,outlying,N=N)

      MUOD_object<-muod(dt, cut_method = c("boxplot", "tangent"))
      MUOD_rst<-evaluation(MUOD_object$outliers$magnitude,outlying,N=N)

      seq0 <- seq_transform(dts = dt, sequence = c("T1", "D1"),  depth_method = "erld", erld_type = "two_sided")
      ST_rst<-evaluation(seq0$outliers$T2,outlying,N=N)

      seq2 <- seq_transform2(dts = dt, dts1 = dt1, sequence = c("SK", "DK"),t_knlg=t_knlg, depth_method = "erld", direction = "upper", outlier_rate = outlier_rate)
      DST_rst<-evaluation(seq2$outliers$DK,outlying,N=N)

      evl_table<-data.frame(rbind(evl_table,c(msplot_rst,TVD_rst,MUOD_rst,ST_rst,DST_rst)))
    }
    names(evl_table)<-paste0(c("msplot_pc","msplot_pf","TVD_pc","TVD_pf","MUOD_pc","MUOD_pf","ST_pc","ST_pf","DST_pc","DST_pf"))
    print(outlier_rate)
    print(data.frame(colMeans(evl_table)))
}
