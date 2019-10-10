
# RustWatch model version from Mika Tode and continue edited by Ying on 10/10/2019
# Shiny application for modelling the evolution of pathogene populations

library(ggplot2)
library(reshape2)
library(Hmisc)
library(xlsx)
library(shinyMatrix)
# library(ggthemes)
#library(D3TableFilter)
library(shiny)
# library(rhandsontable)
# #shiny::runGitHub("rhandsontable", "jrowen",
#  #                subdir = "inst/examples/rhandsontable_output")
library(knitr)
`%then%` <- shiny:::`%OR%`
#load( file = 'r14ken_genind.RData')
#r14ken_inform<-informloci(r14ken_genind)
#r14ken_missing<-missingno(r14ken_inform,cutoff=0)

# This function extracts numbers from a string (e.g. "1,2,3" to 1 2 3)

numextractall <- function(string){ # http://stackoverflow.com/questions/19252663/extracting-decimal-numbers-from-a-string
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)), use.names=FALSE)
}

# Assign shinyServer as an object (here function), in order to use it in the function shinyApp())
shinyServer <- function(input, output, session) {
  
  # This function creates a matrix containing all combinations of 0s and 1s with 2^n rows and n columns
  bincode <- function(n)do.call(expand.grid,rep(list(0:1),n))
  
  #################### Fitness cost vector ############################
  # Creating output for the fitness cost UI, this way it is possible to interact with the fitness cost panel
  # (reactive to n), otherwise it would be static
  
  output$fitnesscostInputUi <- renderUI({
    # Changing in n has to imply a change in the lenght of the fitness cost vector
    n <- input$loci.n
    # Creating a sample vector with length n
    h <- round(runif(n), 2)
    # Code for the UI display
    textInput("fitness.cost","Input fitness costs", value = paste0(h, collapse  = " ") )
  }
  )
  
  
  # Creating the output in the main panel of the fitness cost UI
  # Counts the number input and compares it to the needed variable length n
  
  output$fitness.vector.count <-renderText({
    n <- input$loci.n
    input.text <- as.numeric(numextractall(input$fitness.cost))
    l <- length(input.text)
    if(n>l){
      text.to.show <- paste0('Your input: ',length(input.text), ' number(s). You still need ',n-length(input.text),' number(s) more.')
    }else if(n<l){
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Remove ', length(input.text)-n,' number(s)')
    }
    else{
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Input seems fine.')
    }
  })
  
  ######################### Aggressiveness vector ########################
  # Same procedure as before, main difference: the length of the vectors differ
  
  output$ai.vector.inputUi <- renderUI({
    n <- input$loci.n
    ah <- rep(0.5,2^n)
    textInput("ai_vector","Input for aggressiveness vector ", value = paste0(ah, collapse  = " ") )
  }
  )
  
  output$ai.vector.count <-renderText({
    n <- input$loci.n
    input.text <- as.numeric(numextractall(input$ai_vector))
    l <- length(input.text)
    if(2^n>l){
      text.to.show <- paste0('Your input: ',length(input.text), ' number(s). You still need ',2^n-length(input.text),' number(s) more.')
    }else if(2^n<l){
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Remove ', length(input.text)-2^n,' number(s)')
    }
    else{
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Input seems fine.')
    }
  })
  
  
  
  ######################### jxR matrix ######################################
  # row and column names have to be specified using the matrixInput function
  # ^have been added, see sjxt matrix (SC)
  
  output$jxr.matrix.inputUi <- renderUI({
    n <- input$loci.n
    m <- input$host.m
    vec <- 1:input$host.m
    vec2 <- 1:input$loci.n
    cv <- paste0("Season ", vec)
    rv <- paste0("Locus ", vec2)
    # Function used from the package "ShinyMatrix", useful for creating interactive matrices
    matrixInput("jxr.matrix", value = matrix(sample(c(1,.5,.25,.3,.4),size=m*n,replace = T),nrow=m,ncol=n, dimnames = list(c(cv), c(rv))),rows = list( names = TRUE), cols = list ( names = TRUE), class = "numeric",
                copy = TRUE, paste = TRUE)
  })
  
  
  ######################### sjxt matrix ################################
  # row and column names have to be specified using the matrixInput function,
  # essentially the same as before
  
  # ^have been added with: (SC)
  # dimnames = list(c(cv), c(rv)) 
  # rows = list( names = TRUE), cols = list ( names = TRUE) 
  # vec <- 1:input$seasons.t
  # vec2 <- 1:input$host.m
  # cv <- paste0("Season ", vec)
  # rv <- paste0("Host ", vec2)
  
  # the rows of this matrix have to add up to one in order to contain fractions
  # (it may be good to implement a warning if the user wants to use negative numbers or numbers greater
  # than 1)
  
  output$sjxt.matrix.inputUi <- renderUI({
    t <- input$seasons.t
    m <- input$host.m
    vec <- 1:input$seasons.t
    vec2 <- 1:input$host.m
    cv <- paste0("Season ", vec)
    rv <- paste0("Cultivar ", vec2)
    matrixInput("sjxt.matrix", value = matrix(sample(1/m , size = t*m, replace = T), dimnames = list(c(cv), c(rv)), nrow = t, ncol = m), rows = list( names = TRUE), cols = list ( names = TRUE), class = "numeric",
                copy = TRUE, paste = TRUE)
  })
  
  
  
  ##################### Partial resistance vector (pj vector)################
  # Same procedure for vectors as before
  
  output$pj.vector.inputUi <- renderUI({
    m <- input$host.m
    ap <- rep(0.5,m)
    textInput("pj_vector","Input for the partial resistance vector ", value = paste0(ap, collapse  = " ") )
  }
  )
  
  #
  
  output$pj.vector.count <-renderText({
    m <- input$host.m
    input.text <- as.numeric(numextractall(input$pj_vector))
    l <- length(input.text)
    if(m>l){
      text.to.show <- paste0('Your input: ',length(input.text), ' number(s). You still need ', m-length(input.text),' number(s) more.')
    }else if(m<l){
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Remove ', length(input.text)- m,' number(s)')
    }
    else{
      text.to.show <- paste0('Your input: ',length(input.text), ' numbers. Input seems fine.')
    }
  })
  
  
  
  
  
  
  
  ######################### ixV1 matrix ########################
  
  
  
  output$ixV <-renderTable({
    
    vec <- 1:input$loci.n
    vec2 <- 1:2^input$loci.n
    # Creating row and column names
    cv <- paste0("Locus ", vec)
    rv <- paste0("Phenotype ", vec2)
    # Using the matrix as a dataframe in order to use it with renderTable
    ixV1 <- as.data.frame(bincode(input$loci.n))
    # Renaming 0 with "A" and 1 with "V" for Avirulent and virulent
    ixV1[ixV1 == 0] <- "A"
    ixV1[ixV1 == 1] <- "V"
    colnames(ixV1) <- cv
    rownames(ixV1) <- rv
    ixV1},
    # Showing the rownames in the table, default is FALSE
    rownames = TRUE)
  
  
  
  ######## fi_initial vector with different distributions ##############
  
  # Using renderTable for a vector for a better output
  # Uniform distribution
  output$fi <- renderTable({
    if(input$distr == "unif"){
      fi1 <-  rep(1/(2^input$loci.n), 2^input$loci.n)
      fi1
    }
    # Exponential distribution
    else if (input$distr == "exp"){
      fi1 <-exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))/sum(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum)))
      fi1}
    # Double exponential distribution
    else if (input$distr == "exp2"){
      fi1 <-exp(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum)))/sum(exp(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))))
      fi1
    }
    # Quadratic exponential distribution
    else if (input$distr == "qexp"){
      fi1 <-exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))^2/sum(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))^2)
      fi1
    }
    # Information about the vector in order to have a good display in the UI
    fi1 <- as.data.frame(fi1)
    colnames(fi1) <- "Distribution"
    vec2 <- 1:2^input$loci.n
    rv <- paste0("Phenotype ", vec2)
    rownames(fi1) <- rv
    fi1
  },
  rownames = TRUE)
  
  
  
  
  
  # Creating reactive values, not necessary to look through
  
  
  
  fname = tempfile()
  values = reactiveValues()
  values = reactiveValues()
  
  data = reactive({
    if (!is.null(input$hot)) {
      DF = hot_to_r(input$hot)
    } else {
      if (is.null(values[["DF"]]))
        DF = data.frame(val = 1:10,  nm =rnorm(10),
                        dt = rnorm(10),
                        stringsAsFactors = F)
      else
        DF = values[["DF"]]
    }
    
    
    values[["DF"]] = DF
    DF
  })
  
  
  # Leftover from the old version, pre specified scenarios can still be uploaded from browse and upload
  # This code may be deleted, but does not influence the current version
  output$inputScenarioFiles <- renderUI({
    file.list <- list.files('inputfiles',pattern = '*.xlsx')
    selectInput('scenario','Pre specified scenarios',file.list,selected = file.list[1])
  })
  
  
  #######################  read inputs from a file ###########################
  input.values.from.file <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    
    # Gives a message to the user while the operating system is executing its work
    withProgress(message = 'Reading Data', value = 0, {
      
      # This input is normally used when the user wants to create own input
      # In this version every input is manually changeable and reactive
      
      if (input$readFileInputs=='Use Shiny input'){
        t <- input$seasons.t
        m <- input$host.m
        n <- input$loci.n
        # Saving the "basic" parameters in a dataframe
        parameters <- data.frame(n=n,m=m,t=t)
        cv_vector <- as.numeric(numextractall(input$fitness.cost))
        names(cv_vector) <- paste0('Loci',1:input$loci.n)
        ai_vector <- as.numeric(numextractall(input$ai_vector));names(ai_vector) <- paste0("Aggressiveness", 1:2^n)
        jxr_matrix <- input$jxr.matrix
        pj_vector <- as.numeric(numextractall(input$pj_vector))
        sjxt_matrix <- input$sjxt.matrix
        ixv_matrix <- bincode(n);ixv_matrix <- as.matrix(ixv_matrix)
        # This vector is calculated via a function given in the documentation
        Ci_vector<-(1-apply(1-t(t(ixv_matrix)*cv_vector),1,FUN=prod))
        # Probably already calculated before, it may be possible to delete this part
        # And write fi_initial_vector <- fi1 or fi_initial_vector <- input$fi1
        if(input$distr == "unif"){
          fi1 <-  rep(1/(2^input$loci.n), 2^input$loci.n)
          fi1
        }
        else if (input$distr == "exp"){
          fi1 <-exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))/sum(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum)))
          fi1}
        else if (input$distr == "exp2"){
          fi1 <-exp(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum)))/sum(exp(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))))
          fi1
        }
        else if (input$distr == "qexp"){
          fi1 <-exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))^2/sum(exp(input$loci.n-apply(bincode(input$loci.n),1,FUN=sum))^2)
          fi1}
        
        fi_initial_vector<- fi1
      }
      
      # this part is created for the upload of excel sheets with a certain structure
      # does not currently work!
      # "row.names should be specified
      # There is a problem with a variable and its rownames,
      # I suspect a matrix from the excel sheet
      # It is a bit weird, as the imported rownames are all set to FAlSE
      # Rownames could be needed somewhere though, check differences between Shiny input
      # and the excel sheets
      
      
      else if(input$readFileInputs=='Upload from PC'){
        #return(NULL)
        inFile <- input$file1
        if(is.null(inFile)){return(NULL)}
        parameters <- read.xlsx(inFile$datapath,'parameters')
        n = unlist(parameters['n'])
        cv_vector <- unlist(read.xlsx(inFile$datapath,'cv_vector',header = T,row.names=F))
        #cv_vector <- unlist(read.xlsx(inFile$datapath,'cv_vector',header = T,row.names=T))
        ai_vector <- unlist(read.xlsx(inFile$datapath,'ai_vector',header = T, row.names=F))
        #ai_vector <- unlist(read.xlsx(inFile$datapath,'ai_vector',header = T))
        jxr_matrix<- read.xlsx(inFile$datapath,'jr_matrix',header = T, row.names = F)
        #jxr_matrix<- read.xlsx(inFile$datapath,'jr_matrix',header = T)
        pj_vector <- unlist(read.xlsx(inFile$datapath,'pj_vector',header = T,row.names=F))
        #pj_vector <- unlist(read.xlsx(inFile$datapath,'pj_vector',header = T,row.names=F))
        sjxt_matrix <- read.xlsx(inFile$datapath,'sjxt_matrix',header = T, row.names = T)
        ixv_matrix <- bincode(n); ixv_matrix <- as.matrix(ixv_matrix)
        #hans fi initial vector
        # The distribution of the fi initial vector has to be deducted from the excel sheet too
        # fi_initial_vector <- read.xlsx(inFile$datapath,'fi_initial_vector',header = F)
        # Just delete the line below and it will work, ONLY IF the excel sheet contains information
        # about this vector and a sheet called "fi_initial_vector"!
        fi_initial_vector<-exp(n-apply(ixv_matrix,1,FUN=sum))/sum(exp(n-apply(ixv_matrix,1,FUN=sum)))
        Ci_vector<-(1-apply(1-t(t(ixv_matrix)*cv_vector),1,FUN=prod))
      }
      
      # Save the parameters in a list (independent from the method chosen above)
      
      input.values <- list(parameters=parameters,ixv_matrix=ixv_matrix,cv_vector=cv_vector,Ai_vector=ai_vector,jxR_matrix=jxr_matrix,Pj_vector=pj_vector,fi_initial_vector=fi_initial_vector,sjxt_matrix=sjxt_matrix,Ci_vector=Ci_vector)
      incProgress(1, detail = paste("Done"))})
    if(!exists("input.values")){return(NULL)}
    return(input.values)
    
  })
  
  # Is not used anymore, can be used in order to show a nice table for the "basic" parameters
  # Some modification in the UI has to be done then
  output$sample.parameter <-  renderTable({
    file.parameters <- input.values.from.file()
    if(is.null(file.parameters)){
      return(NULL)
    }
    parameters <- data.frame(file.parameters[[1]])
    
    return(parameters)
    
  })
  
  # Same here, shows selected parameters (jxr matrix, sjxt matrix, etc.)
  # Unnessecary as it is already shown in the input
  # Maybe later summarised in the output window?
  
  output$fileParameter <- DT::renderDataTable({
    parameter.selected <- input$sel.file.parameter
    parameters.from.file <- input.values.from.file()
    parameter.value <- data.frame(parameters.from.file[[parameter.selected]])
    return(parameter.value)
  })
  
  
  # Model computations, do not change
  
  computations1 <- reactive({
    
    # unlist all parameters from the respective data
    # The mistake for "rownames must be specified" could be here, check what "unlist" does
    parameters.from.file <- input.values.from.file()
    jxR_matrix<-as.matrix(parameters.from.file[['jxR_matrix']])
    ixV_matrix<-as.matrix(parameters.from.file[['ixv_matrix']])
    Pj_vector <- parameters.from.file[['Pj_vector']]
    Ci_vector<- parameters.from.file[['Ci_vector']]
    Ai_vector<- parameters.from.file[['Ai_vector']]
    n <- unlist(parameters.from.file[['parameters']]['n'])
    m <- unlist(parameters.from.file[['parameters']]['m'])
    ujxi_matrix<-matrix(NA,nrow=m,ncol=2^n)
    pp<-rep(NA,n)
    # Creating the "information" matrix uxji based on several inputs
    for (j in 1:m){
      for (i in 1:2^n){
        for (k in 1:n) {
          pp[k]<-min(1,(1-jxR_matrix)[j,k]+ixV_matrix[i,k])
        }
        ujxi_matrix[j,i]<-min(1,pp)%*%((1-Pj_vector[j])*(1-Ci_vector[i])*Ai_vector[i])
      }
    }
    return(ujxi_matrix)
  })
  
  
  ############################# ujxi matrix ########################
  
  output$ujxi <-renderTable({
    n <- input$loci.n
    m <- input$host.m
    v1 <-  1:m
    # Rownames
    cv <- paste0("Cultivar ", v1)
    # Using the reactive computations
    ujxi <- computations1()
    # Creating the colnames (AAA, AAV, etc.)
    nf<-bincode(n)
    nf[nf==0]<-'A'
    nf[nf==1]<-'V'
    lv <-  numeric(0)
    for (i in 1:2^n){
      lv <- c(lv, paste(nf[i,],collapse = ""))
    }
    ujxi <- as.data.frame(ujxi)
    rownames(ujxi) <- cv
    colnames(ujxi) <- lv
    # This last line has to be included in order to have an output
    ujxi
  },
  rownames = TRUE)
  
  
  # Again: some computations, rownames may be used somewhere in the code and create the error
  ##################### computation 2 #############
  computation2 <- reactive({
    parameters.from.file <- input.values.from.file()
    fi_initial_vector<- unlist(parameters.from.file[['fi_initial_vector']])
    sjxt_matrix <- as.matrix(parameters.from.file[['sjxt_matrix']])
    tmax <- unlist(parameters.from.file[['parameters']]['t'])
    n <- unlist(parameters.from.file[['parameters']]['n'])
    ixV_matrix<-as.matrix(bincode(n))
    # computation1.results <- computations1()
    ujxi_matrix<-  computations1()
    wjt_vector<-fi_initial_vector%*%t(ujxi_matrix)
    wjsj<-wjt_vector*sjxt_matrix[1,]
    #sum(wjsj)
    f_prime_<-data.frame(season_1=fi_initial_vector*apply(ujxi_matrix*sjxt_matrix[1,],2,sum)/sum(wjsj))
    w_prime_<-c(NA,sum(wjsj))
    validate(
      need(!is.null(parameters.from.file), "File Input is required. Choose file and upload for results")  %then%
        need(!is.null( ujxi_matrix), "File Input is required. Choose file and upload for results")
    )
    for (t in 2:tmax){
      wjt_vector<-f_prime_[,(t-1)]%*%t(ujxi_matrix)
      wjsj<-wjt_vector*sjxt_matrix[t,]
      f_prime_<-cbind(f_prime_,as.numeric(f_prime_[,(t-1)])*apply(ujxi_matrix*sjxt_matrix[t,],2,sum)/sum(wjsj))
      names(f_prime_)[t]<-paste0("season_",t)
      w_prime_<-c(w_prime_,sum(wjsj))
    }
    namesf_prime_<-ixV_matrix
    namesf_prime_[namesf_prime_==0]<-'A'
    namesf_prime_[namesf_prime_==1]<-'V'
    namesf_prime_vector<-c()
    for (na in 1:dim(ixV_matrix)[1]){
      ppp<-namesf_prime_[na,1]
      for(xy in 2:dim(ixV_matrix)[2]){
        ppp<-paste0(ppp,namesf_prime_[na,xy])
      }
      namesf_prime_vector<-c(namesf_prime_vector,ppp)
    }
    
    f_prime_<-t(cbind(season_0=fi_initial_vector,f_prime_))
    colnames(f_prime_)<-namesf_prime_vector
    
    
    virulence <- f_prime_%*%ixV_matrix
    
    colnames(virulence)<-paste0("V",1:n)
    
    ################################################################################################################
    # Binding of matrices to facilitate outputs
    f_prime_w_prime_<-rbind(w_prime_,t(f_prime_))
    
    
    names(w_prime_)<-rownames(f_prime_)
    
    
    w_prime_melt<-melt(w_prime_)
    w_prime_melt$season<-as.numeric(factor(rownames(w_prime_melt),levels=paste0("season_",0:tmax),labels=0:tmax))-1
    
    virulence_melt<-melt(virulence)
    virulence_melt$season<-as.numeric(factor(virulence_melt[,1],levels=paste0("season_",0:tmax),labels=1:length(levels(virulence_melt[,1]))))-1
    
    names(virulence_melt)[2]<-"virulence"
    virulence_w_prime_melt<-merge(virulence_melt,w_prime_melt, by="season")
    names(virulence_w_prime_melt)[4:5] <- c("relfre", "meanfitness")
    
    #rownames(f_prime_)[1:dim(ixV_matrix)[1]]<-namesf_prime_vector
    #f_prime_<-as.data.frame(cbind(f_prime_,Loci=namesf_prime_vector))
    
    f_prime_melt<-melt(f_prime_)
    f_prime_melt$season<-as.numeric(factor(f_prime_melt[,1],levels=paste0("season_",0:tmax),labels=1:length(levels(f_prime_melt[,1]))))-1
    names(f_prime_melt)[1:2]<-c("seas","Loci")
    
    w_prime_melt$seas <- row.names(w_prime_melt)
    f_prime_w_prime_melt<-merge(f_prime_melt,w_prime_melt, by=c("season","seas"))
    names(f_prime_w_prime_melt)[4:5] <- c("relfre", "meanfitness")
    
    sjxt_matrix_melt<-melt(sjxt_matrix)
    names(sjxt_matrix_melt)<-c("season","host_genotype","area_fraction")
    sjxt_matrix_melt$seas<-factor(sjxt_matrix_melt[,1],labels=paste0("season_",1:tmax),levels=1:length(sjxt_matrix[,1]))
    sjxt_matrix_w_prime_melt<-merge(sjxt_matrix_melt,w_prime_melt, by=c("season","seas"))
    names(sjxt_matrix_w_prime_melt)[5]<-"meanfitness"
    sjxt_matrix_w_prime_melt$host_genotype<-as.factor(sjxt_matrix_w_prime_melt$host_genotype)
    
    
    
    computation2.results <- list(w_prime_=w_prime_,sjxt_matrix=sjxt_matrix,virulence_w_prime_melt=virulence_w_prime_melt,f_prime_w_prime_melt=f_prime_w_prime_melt,ujxi_matrix=ujxi_matrix,sjxt_matrix_w_prime_melt=sjxt_matrix_w_prime_melt)
    return(computation2.results)
  })
  
  ####################################### Plots #################################
  
  ##  # mean relative fitness of pathogen population & relative host genotype areas over time
  output$plot1 <- renderPlot({
    #   validate(
    #     need(!is.null(input$file1), "File Input is required. Choose file and upload for results")
    #   )
    #  if (is.null(input$file1))
    #   return(NULL)
    
    withProgress(message = 'Generating Plot ...', value = 0, {
      # parameters are saved
      parameters.from.file <- input.values.from.file()
      computation2.results <- computation2()
      tmax <- unlist(parameters.from.file[['parameters']]['t'])
      m <- unlist(parameters.from.file[['parameters']]['m'])
      n <- unlist(parameters.from.file[['parameters']]['n'])
      sjxt_matrix <- computation2.results[['sjxt_matrix']]
      w_prime_<- computation2.results[['w_prime_']]
      ujxi_matrix <- computation2.results[['ujxi_matrix']]
      # Maybe replacing 0 and 1 here too?
      ixv_matrix<-bincode(n)
      fi_initial_vector<-exp(n-apply(ixv_matrix,1,FUN=sum))/sum(exp(n-apply(ixv_matrix,1,FUN=sum)))
      sjxt_matrix_w_prime_melt <- computation2.results[['sjxt_matrix_w_prime_melt']]
      
      # Actual plot
      # Line size was reduced, the overall layout was improved, help texts should be
      # further specified in the UI
      p <-  ggplot(sjxt_matrix_w_prime_melt) + geom_line(aes(x=season, y=area_fraction, colour=host_genotype),size=rel(0.9))+
        geom_line(aes(x=season, y=meanfitness),color="black",size=rel(0.7),linetype=2)  +
        labs(y="Host genotype area fraction resp. mean fitness",x="Season")+
        scale_x_continuous(limits=c(1,(tmax)),breaks=1:(tmax), labels=1:(tmax))+
        scale_y_continuous(limits=c(0,1))+
        theme_classic()
      incProgress(1, detail = paste("Done"))})
    return(p)
  })
  
  ######################### Second plot ################################
  
  generatePlot2 <- reactive({
    plotnumber <- input$plot2.no
    computation2.results <- computation2()
    virulence_w_prime_melt <- computation2.results[['virulence_w_prime_melt']]
    f_prime_w_prime_melt <- computation2.results[['f_prime_w_prime_melt']]
    parameters.from.file <- input.values.from.file()
    tmax <- unlist(parameters.from.file[['parameters']]['t'])
    p <- ggplot(virulence_w_prime_melt) + geom_line(aes(x=season, y=relfre, colour=virulence),size=rel(0.9))+ geom_point(aes(x=season, y=relfre, colour=virulence),size=rel(1.2)) +
      geom_line(aes(x=season, y=meanfitness),color="black",size=rel(0.7),linetype=2)  +
      labs(y="Relative frequency virulence properties",x="Season")+
      scale_x_continuous(limits=c(0,(tmax)),breaks=0:(tmax), labels=0:(tmax))+
      theme_classic()
    return(p)
    
  })
  
  
  ######################## Third plot #############################
  generatePlot3 <- reactive({
    computation2.results <- computation2()
    virulence_w_prime_melt <- computation2.results[['virulence_w_prime_melt']]
    f_prime_w_prime_melt <- computation2.results[['f_prime_w_prime_melt']]
    parameters.from.file <- input.values.from.file()
    tmax <- unlist(parameters.from.file[['parameters']]['t'])
    p <-  ggplot(f_prime_w_prime_melt) + geom_line(aes(x=season, y=relfre, colour=Loci),size=rel(0.9))+
      geom_line(aes(x=season, y=meanfitness),color="black",size=rel(0.7),linetype=2)  +
      labs(y="Relative frequency resp. mean fitness",x="Season")+
      scale_x_continuous(limits=c(0,(tmax)),breaks=0:(tmax), labels=0:(tmax))+
      scale_y_continuous(limits=c(0,1))+
      theme_classic()
    
    return(p)
    
  })
  
  # Creating reactive outputs
  
  output$plot2 <- renderPlot({
    withProgress(message = 'Generating Plot', value = 0, {
      
      p <- generatePlot2()
      
      incProgress(1, detail = paste("Done"))})
    return(p)
  })
  
  output$plot3 <- renderPlot({
    withProgress(message = 'Generating Plot', value = 0, {
      
      p <- generatePlot3()
      
      incProgress(1, detail = paste("Done"))})
    return(p)
  })
  
  # Creating a downloadable Scenario file
  # Does not currently work!
  # Probably an error addressing the names of the used data files
  
  output$report = downloadHandler(
    
    
    filename = 'report/myreport.pdf',
    
    content = function(file) {
      out = knit2pdf('report/input.Rnw', clean = F)
      file.rename(out, file) # move pdf to file for downloading
    },
    
    contentType = 'application/pdf'
  )
  
  # Not useful, as this creates random parameters and saves them
  # Chosen vectors and matrices should be included here
  # Basically use the list created before the computations should help
  ## making a sample file with parameters.
  make.sample.parameter.file <- reactive({
    t <- input$seasons.t
    m <- input$host.m
    n <- input$loci.n
    parameters <- data.frame(n=n,m=m,t=t)
    #ixv_matrix <- bincode(n)
    cv_vector <- runif(n);names(cv_vector) <- paste0('LOCI',1:n)
    cv_vector <- data.frame(cv_vector)
    ai_vector <- data.frame(aggressiveness=rep(0.5,2^n))
    jr_matrix <- matrix(sample(c(1,.5,.25,.3,.4),size=m*n,replace = T),nrow=m,ncol=n)
    pj_vector <- runif(m);pj_vector <- data.frame(pj_vector)
    row.names(pj_vector) <- paste0('Host genotype ',1:nrow(pj_vector))
    #fi_initial_vector<-exp(n-apply(ixV_matrix,1,FUN=sum))/sum(exp(n-apply(ixV_matrix,1,FUN=sum)))
    sjxt_matrix <- matrix( sample(c(.1,.5,.2),size=t*m,replace = T),nrow=t,ncol=m)
    file.to.save <- paste0('report/sampleParameters',n,'_',m,'_',t,'.xlsx')
    write.xlsx(parameters,file.to.save,sheetName='parameters',col.names = T,row.names = F,append = F)
    write.xlsx(cv_vector,file.to.save,sheetName='cv_vector',col.names = T,row.names = F,append = T)
    write.xlsx(ai_vector,file.to.save,sheetName='ai_vector',col.names = T,row.names = F,append = T)
    write.xlsx(jr_matrix,file.to.save,sheetName='jr_matrix',col.names = T,row.names = F,append = T)
    write.xlsx(pj_vector,file.to.save,sheetName='pj_vector',col.names = T,row.names = F,append = T)
    write.xlsx(sjxt_matrix,file.to.save,sheetName='sjxt_matrix',col.names = T,row.names = F,append = T)
    
  })
  
  
  output$sample.parameter.file = downloadHandler(
    
    filename = function(){
      t <- input$seasons.t
      m <- input$host.m
      n <- input$loci.n
      filename=paste0('report/sampleParameters',n,'_',m,'_',t,'.xlsx')
      return(filename)
    },
    content <- function(file) {
      t <- input$seasons.t
      m <- input$host.m
      n <- input$loci.n
      make.sample.parameter.file()
      filename=paste0('report/sampleParameters',n,'_',m,'_',t,'.xlsx')
      file.copy(filename, file)
      #      file.rename(file.to.save,file)
    },
    contentType = 'application/xlsx'
  )
  
  
}
# end shiny server

# This has to be run AFTER both shinyServer() and shinyUI were executed
# Creates the actual shiny application

#shinyApp(ui = shinyUI, server = shinyServer)

