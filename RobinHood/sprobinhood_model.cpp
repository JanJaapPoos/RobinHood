
  #include <TMB.hpp>
  using namespace density;
  template <class Type>
  Type objective_function<Type>::operator () (){
    //------
    // DATA
    //------
    // data matrices
    DATA_MATRIX(Imat);
    DATA_MATRIX(Cmat);
    DATA_MATRIX(Fmat);
    DATA_IMATRIX(nIndices);
    DATA_IMATRIX(Ipresent);
    DATA_MATRIX(corrmatrix);
    DATA_MATRIX(sdrwdr); // sds for random walks of Fvalues transposed to fractions like U coming from data rich
    //------------
    // PARAMETERS
    //------------
      // put in parameter for corrmult somewhere
    PARAMETER_VECTOR(logm);
    PARAMETER_VECTOR(logK);
    PARAMETER_VECTOR(logq);
    PARAMETER_VECTOR(logsdi);
    PARAMETER_VECTOR(logsdc);
    PARAMETER_MATRIX(logsdrw); // log sds estimated for random walks of data poor
    PARAMETER_VECTOR(logitalpha);
    PARAMETER_MATRIX(logitU);
    // exploitation rate (U of stock assessments)
    
    
    //---------------------------
    // PRELIMINARY CALCULATIONS
    //---------------------------
    // NUMBER OF OBS IN CATCHES AND F
    int n = Imat.rows();
    // NUMBER OF STOCKS IN UNSTRUCTURED CATCHES
    int mc = Cmat.cols();
    // NUMBER OF INDICES
    int mi = nIndices.sum();
    // NUMBER OF STOCKS IN EXTERNAL Fs
    int mf = Fmat.cols();
    // Total number of stocks for which we have F/U info
    int mt = mc + mf;
    
    
    // rework parms
    std::cout << mi << std::endl;
  
    // put in parameter for corrmult somewhere
    vector<Type> m = exp(logm);
    vector<Type> K = exp(logK);
    vector<Type> q = exp(logq);
    vector<Type> sdi = exp(logsdi);
    vector<Type> sdc = exp(logsdc);
    vector<Type> alpha(mc);
     for(int j = 0; j < mc; j++){
      alpha(j) = 1.0 / (1.0 + exp(-logitalpha(j))); // alpha on (0,1)
    }
    
    matrix<Type> sdrw = logsdrw.array().exp();
    
    
    // container for fitted index and catches, and for all sd rws (estimated for data poor and taken for DR
    matrix<Type> Ihat(n, mi); // we have a prediction for the actual I with sigma that construct I_observations
    matrix<Type> Chat(n, mc); // likewise, but this one also involves the rules for the random walk
    matrix<Type> sdrwall(mt,1); // these sigmas will later be used in the multivariate normal distribution with a correlation parameter to create the jumps in time for U1

    
    
    // fill container all sd rw with default sd for U1
    for(int j = 0; j < mc; j++)
      sdrwall(j,0) = sdrw(j);
    for(int j = 0; j < mf; j++)
      sdrwall(mc+j,0) = sdrwdr(j);
    
    //-----------
    // PROCEDURE
    //-----------
    Type nll = 0.0; // initialize
    matrix<Type> B(n, mc);
    // start at K * alpha

    matrix<Type> U = 1.0 / (1.0 + (-logitU).array().exp());
   
    // get better eigen fill //what is ment here?
    for(int j = 0; j < mc; j++){
      B(0,j) = alpha(j) * K(j);  //B0 relies on chosen alpha
    }
    for(int i = 1; i < n; i++){ 
      for(int j = 0; j < mc; j++){
        //B(i,j) = B(i-1,j) +   4.*  m(j)  * (B(i-1,j)/K(j))  - 4.*  m(j)  * pow( ( B(i-1,j) / K(j)),2) - (U(i-1,j) * B(i-1,j));
        B(i,j) = abs(B(i-1,j) +   4.*  m(j)  * (B(i-1,j)/K(j))  - 4.*  m(j)  * pow( ( B(i-1,j) / K(j)),2) - (U(i-1,j) * B(i-1,j)));
      }
    }
    
    matrix<Type> Sigma = (sdrwall * (sdrwall.transpose())).array() * corrmatrix.array(); //magical spell that i still have to figure out
    // if i would add corrmult where would i put it? -> matrix<Type> Sigma = (sdrwall * (sdrwall.transpose())).array() * variable_corrmatrix.array();
      //create extra line of code that calculates new corrmatrix from original one and optimize that one-> 
      // matrix<type>variable_corrmatrix.array()= corrmult* corrmatrix.array()
      //for (int(variable_corrmatrix.array() blablabla function that transforms the diagonal to 1
      
    //the "matrix<type>" line takes the default sigmas and trasformes them with a covariance matrix
    //these are calculated from the correlation matrix 
    //here is the sigma of the process, but where is it in the R  
    for(int i = 1; i < n; i++){

      nll += MVNORM(Sigma)(vector<Type>(logitU.row(i) - logitU.row(i-1))); // note + as return negative log density
      //if I am right, this is the process part of the random walk for the Fishing mortalities
      //MVNORM is similar to dnorm
      
      
      //std::cout << std::endl << std::endl;
      //std::cout  << "nll " << nll  << std::endl << std::endl;

    }
    // CM
    // observations component for data rich F to make the logU's follow the observations
    for(int i = 0; i < n; i++){ 
      for(int j = 0; j < mf; j++){
   	    nll -= dnorm(log(Fmat(i,j)),log(U(i,(mc + j))), Type(1.0e-3), true);
      }
    }
    //so this is then the observation part, where we fit the estimates for F to the data we have
    //what is F mat: only data from data rich here
    //why is it here suddenly only logU instead of logitU?
    
    
   //std::cout  << "nll after Fmat " << nll ;
   //std::cout << std::endl << std::endl;
        
    // B/BMSY
    matrix<Type> BoBMSY(n,mc);
    matrix<Type> FoFMSY(n,mc);

    for(int i = 0; i < n ; i++){ 
      for(int j = 0; j < mc; j++){
        BoBMSY(i,j)= B(i,j)/(K(j) *0.5);
	      FoFMSY(i,j)= U(i,j)/ (m(j)/(K(j) * 0.5));
      }
    }
    // LIKELIHOOD
    // observations component indices data poors (mc = num stocks, n = num years)
    int idx = 0;
    for(int j = 0; j < mc; j++){
      for(int k = 0; k < nIndices(j); k++){
	      std::cout << idx << std::endl;
	      for(int i = 0; i < n; i++){ 
	        Ihat(i,idx) = q(idx) * B(i,j);
	        // check if index observation is present
	        if(Ipresent(i,idx) == 1){
	            nll -= dnorm(log(Imat(i,idx)), log(Ihat(i,idx)) - pow(sdi(idx), 2.0) / 2.0, sdi(idx), true); 
	        }		     
	      }
	      idx += 1;
      }
    }
    std::cout << std::endl << std::endl;
    std::cout  << "nll after Ihat " << nll ;
    std::cout << std::endl << std::endl;
    
    
    //observation component catches data poors
    for(int i = 0; i < n-1; i++){ 
      for(int j = 0; j < mc; j++){
	      Chat(i,j) = U(i,j) * B(i,j);
   	    nll -= dnorm(log(Cmat(i,j)),log(Chat(i,j)), sdc(j), true);
      }
    }

    std::cout << std::endl << std::endl;
    std::cout  << "nll after Chat " << nll ;
    std::cout << std::endl << std::endl;
    
    ADREPORT(B);
    ADREPORT(U); 
    ADREPORT(BoBMSY);
    ADREPORT(FoFMSY);
    ADREPORT(Ihat);
    ADREPORT(Chat);
    ADREPORT(Sigma);
    return nll;
  }
  
  //it is a giant but the assumption is that the lowest nll will eject the best parameters for the model. 
  //assumptions include that errors are lognomrally distributed
  
  //also wonder why one would not only find the parameters of the part of the model that is most important first via nll 
  //and then from there on work out the rest of the parameters (also when overall nll is higher with this method)

