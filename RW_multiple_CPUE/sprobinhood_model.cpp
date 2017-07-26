
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
    DATA_MATRIX(sdrwdr); // sds for random walks of U/F coming from data rich
    //------------
    // PARAMETERS
    //------------
    PARAMETER_VECTOR(logm);
    PARAMETER_VECTOR(logK);
    PARAMETER_VECTOR(logq);
    PARAMETER_VECTOR(logsdi);
    PARAMETER_VECTOR(logsdc);
    PARAMETER_MATRIX(logsdrw); // log sds estimated for random walks of data poor
    PARAMETER_VECTOR(logitalpha);
      // exploitation rate (U of stock assessments)
    PARAMETER_MATRIX(logitU);
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
    matrix<Type> Ihat(n, mi);
    matrix<Type> Chat(n, mc);
    matrix<Type> sdrwall(mt,1);

    // fill container all sd rw
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

    matrix<Type> U = 2.0 / (1.0 + (-logitU).array().exp());
   
    // get better eigen fill
    for(int j = 0; j < mc; j++){
      B(0,j) = alpha(j) * K(j);
    }
    for(int i = 1; i < n; i++){ 
      for(int j = 0; j < mc; j++){
        B(i,j) = abs(B(i-1,j) +   4.*  m(j)  * (B(i-1,j)/K(j))  - 4.*  m(j)  * pow( ( B(i-1,j) / K(j)),2) - (U(i-1,j) * B(i-1,j)));
      }
    }
    matrix<Type> Sigma = (sdrwall * (sdrwall.transpose())).array() * corrmatrix.array();
        
    for(int i = 1; i < n; i++){

      nll += MVNORM(Sigma)(vector<Type>(logitU.row(i) - logitU.row(i-1))); // note + as return negative log density

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
    //observation component catches data poors
    for(int i = 0; i < n-1; i++){ 
      for(int j = 0; j < mc; j++){
	Chat(i,j) = U(i,j) * B(i,j);
   	nll -= dnorm(log(Cmat(i,j)),log(Chat(i,j)), sdc(j), true);
      }
    }

    std::cout << std::endl << std::endl;
    std::cout  << "nll " << nll ;
    std::cout << std::endl << std::endl;
    
    ADREPORT(B);
    //ADREPORT(U);
    ADREPORT(BoBMSY);
    ADREPORT(FoFMSY);
    ADREPORT(Ihat);
    ADREPORT(Sigma);
    return nll;
  }

