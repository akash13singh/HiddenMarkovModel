
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;


public class BaumWelch {

	/*
	 * TODO use setters and getters for these variables.
	 */
	private double[][] transitionMatrix = null;
	private double[][] emmissionMatrix = null;
	private double[][] initialStateDistribution = null;
	private int[] emissionSequence = null;
	//c: scale factor
	private double[] c = null;
	
	// an N*T matrix used in alpha/forward pass
	private double[][] beta = null;
	private double prob_model = 0 ;
	private double[][][] digamma = null;
	// an N*T matrix used in alpha/forward pass
	private double[][] alpha = null;

	private double[][] gamma = null;

	
	// number of states
	private int N;
	
	
	private int K ;
	
	// length of emmission sequence;
	private int T;
	
	/*
	 * reads input file through pipe
	 */
	private  String[] readInput() throws IOException{
		String[] input = new String[4];
		InputStreamReader inputStream = new InputStreamReader(System.in);
		BufferedReader bufferedReader = new BufferedReader(inputStream);
		int i =0;
		while(i<4){
			input[i] = bufferedReader.readLine();
			i++;
		}
		return input;
	}
	
	/*
	 * code to convert a string of numbers separated by space into 2d arrays. 
	 * array dimesnions are the first two numbers in the string
	 */
	private double[][] convertStringTo2DArray(String str){
		String[] arr = str.split(" ");
		int rows = Integer.parseInt(arr[0]);
		int cols = Integer.parseInt(arr[1]);
		int counter = 2;
		double[][] doubleArray= new double[rows][cols];
		for(int i=0;i<rows;i++){
			for(int j=0;j<cols;j++){
				doubleArray[i][j] = Double.parseDouble(arr[counter]);
				counter++;
			}
		}
		return doubleArray;
	}
	
	private int[] convertStringTo1DArray(String str){
		String[] stringArray = str.split(" ");
		int num = Integer.parseInt(stringArray[0]);
		int counter = 1;
		int[] array = new int[num ];
		for(int i=0;i<num ;i++){
			array[i] = Integer.parseInt(stringArray[counter]);
			counter++;
		}
		return array;
	}
	
	/*
	 * code to multiply two matrices
	 */
	private double[][] matrixMultiplication(double[][] a, double[][] b){ 
	  int rows = a.length;
	  int loop= a[0].length;
	  int cols = b[0].length;
	  double[][] c = new double[rows][cols];
	  
	  // uncomment the below lines for debugging
	  /*
	     System.out.println(Arrays.deepToString(a));
	     System.out.println(Arrays.deepToString(b))	;
	  */
	  
	  for(int i=0;i<rows;i++){
		  for(int j=0;j<loop;j++){
			  for(int k =0;k<cols;k++){
				  c[i][k] += a[i][j]*b[j][k];
				  //c[i][k] = new BigDecimal(c[i][k]).setScale(2, RoundingMode.HALF_UP).doubleValue();
			  }
			 
		  }
	  }
	  
	  return c;
	}
	
	private void alphaPass(){
		
	    c = new double[T];
		alpha = new double[N][T];
		
	    //calculate alpha 0;
		c[0]=0 ;
		for(int i = 0; i<N; i++){
			alpha[i][0]= initialStateDistribution[0][i] * emmissionMatrix[i][emissionSequence[0]];
			c[0] = c[0] + alpha[i][0];
		}		
		c[0] = 1/c[0];	
		
		//scale alpha[i][0]
		for(int i = 0; i<N; i++){
			alpha[i][0]= c[0] * alpha[i][0]; 
		}
		
		//compute alpha[i][t]
		for(int t =1; t<T; t++){
			c[t]=0;
			for(int i = 0; i<N; i++){
				alpha[i][t]=0;
				for(int j=0; j<N;j++){
					alpha[i][t] = alpha[i][t]+  alpha[j][t-1] * transitionMatrix[j][i];
				}
				alpha[i][t] = alpha[i][t] * emmissionMatrix[i][emissionSequence[t]];
				c[t]= c[t]+alpha[i][t];
			}
			//System.out.println(Arrays.deepToString(alpha));
			//scale alpha[i][t]
			c[t] = 1/c[t];
			for(int i = 0; i<N; i++){
				alpha[i][t] = c[t]*alpha[i][t];	
			}
		}
	/*	System.out.println(Arrays.deepToString(alpha));
		System.out.println("============");*/
		double prob = 0 ;
		for (int i =0 ; i< N ; i++)
		{
			prob += alpha[i][T-1];
		}
		prob_model = prob ;
	     prob = 1 ;
		for(int t = 0; t<T; t++){
			prob = prob*c[t];
		}
		prob_model =  1/prob;
		//System.out.println(prob_model);
	}
	
	private void betaPass(){	
	   // c = new double[T];
		beta = new double[N][T];
		
	    //calculate beta T-1;
		//c[T-1] =0 ;
		for(int i = 0; i<N; i++){
			beta[i][T-1]= c[T-1] ; 
			//c[T-1] = c[T-1] + beta[i][T-1];
		}		
		//c[T-1] = 1/c[T-1];	
	
		
		//compute beta[i][t]
		for(int t =T-2; t>=0; t--){
			for(int i = 0; i<N; i++){
				beta[i][t]=0;
				for(int j=0; j<N;j++){
					beta[i][t] +=  (beta[j][t+1] * transitionMatrix[i][j] * emmissionMatrix[j][emissionSequence[t+1]]);
				}
				beta[i][t] *= c[t];
			}
		}
	}
	
	private void di_gamma_func(){	
		
		digamma = new double[N][T][N];
		
		for(int t =0; t<T-1; t++){
			for(int i = 0; i<N; i++){
				for(int j=0; j<N;j++){	
				/*	double ii = alpha[t][i] ;
					 ii =  transitionMatrix[i][j]  ;
					 ii = alpha[t][i] ;
					 ii = alpha[t][i] ;
					 ii = emmissionMatrix[j][emissionSequence[t+1]] ;
					 ii =beta[t+1][j];*/
					digamma[i][t][j]= (alpha[i][t]/prob_model) * transitionMatrix[i][j] * emmissionMatrix[j][emissionSequence[t+1]]* beta[j][t+1];
					//System.out.println(digamma[i][t][j]);
				}
			}
		}	
		
		gamma = new double[N][T] ;
		for(int t =0; t<T-1; t++){
			for(int i = 0; i<N; i++){
				gamma[i][t] = 0 ;
				for(int j=0; j<N;j++){	
					gamma[i][t] +=digamma[i][t][j];
				}
			}
		}	
	}
	
	
	
	

	private void di_gamma_func2(){
		
		digamma = new double[N][T][N];
		gamma = new double[N][T] ;
		for(int t =0; t<T-1; t++){
			double denom = 0 ;
			for(int i = 0; i<N; i++){
				for(int j=0; j<N;j++){ 
					denom += alpha[i][t] * transitionMatrix[i][j] * emmissionMatrix[j][emissionSequence[t+1]]* beta[j][t+1];		
				}
			}
			for(int i = 0; i<N; i++){
				gamma[i][t] = 0 ;
				for(int j=0; j<N;j++){ 
					
					digamma[i][t][j]= alpha[i][t] * transitionMatrix[i][j] * emmissionMatrix[j][emissionSequence[t+1]]* beta[j][t+1]/denom;
					gamma[i][t] +=digamma[i][t][j];
				}
			}
		}
		
		
				double denom = 0;
				
				for ( int i =0 ; i < N ; i++){
					
					denom+= alpha[i][T-1];
				}	
	          for ( int i =0 ; i < N ; i++){
	        	  gamma[i][T-1] = alpha[i][T-1]/denom;
				}	
	}
	
	
	private void baum_welch(){
		transitionMatrix = new double[N][N];
		emmissionMatrix = new double[N][K];
		initialStateDistribution = new double[1][N];
		
		for(int i = 0; i<N; i++){ 
			initialStateDistribution[0][i] = gamma[i][0];
			
		}
		
		for(int i = 0; i<N; i++){ 
			for(int j = 0; j<N; j++){ 
				double num_sum = 0 ;
				double den_sum = 0 ;
				for(int t =0; t<T-1; t++){
					num_sum += digamma[i][t][j];
					den_sum += gamma[i][t];
				}
				//System.out.println(den_sum);			
				transitionMatrix[i][j]=num_sum/den_sum;
			}
		}
		
	
		for(int j = 0; j<N; j++){ 
			for(int k = 0; k<K; k++){ 
				double num_sum = 0 ;
				double den_sum = 0 ;
				for(int t =0; t<T; t++){
					if (emissionSequence[t] == k) {
					    num_sum += gamma[j][t];
					}
					den_sum += gamma[j][t];
				}
				emmissionMatrix[j][k]=num_sum/den_sum;
			}
		}
		
		
		
/*		for(int i = 0; i<N; i++){ 
			double row_sum = 0;
			for(int j = 0; j<N; j++){ 
				row_sum += transitionMatrix[i][j];
			}
			for(int j = 0; j<N; j++){ 
				transitionMatrix[i][j]/=row_sum;
			}
			
		}
		
		for(int i = 0; i<N; i++){ 
			double row_sum = 0;
			for(int k = 0; k<K; k++){ 
				row_sum += emmissionMatrix[i][k];
			}
			for(int k = 0; k<K; k++){ 
				emmissionMatrix[i][k]/=row_sum;
			}
			
		}*/
		
/*		System.out.println("====================================!!!!!!!!");
		System.out.println(Arrays.deepToString(A));
		System.out.println(Arrays.deepToString(B));
		System.out.println("====================================!!!!!!!!");*/

		
	}
	
	
	
	private void calculateEmissionProbability(String[] hmmParams){
		double prob =1;
		transitionMatrix = convertStringTo2DArray(hmmParams[0]);
		emmissionMatrix = convertStringTo2DArray(hmmParams[1]);
		initialStateDistribution = convertStringTo2DArray(hmmParams[2]);
		emissionSequence = convertStringTo1DArray(hmmParams[3]);
		
//		uncomment for debugging
//		System.out.println(Arrays.deepToString(transitionMatrix));
//		System.out.println(Arrays.deepToString(emmissionMatrix));
//		System.out.println(Arrays.deepToString(initialStateDistribution));
//		System.out.println(Arrays.toString(emissionSequence));
		
		N = transitionMatrix.length;
		K= emmissionMatrix[0].length;
		T= emissionSequence.length;
		for (int i =0 ; i<15; i++)
		{
		alphaPass();		
		betaPass();
		di_gamma_func2();
	    System.out.println(prob_model);
		baum_welch();
		System.out.println(Arrays.deepToString(transitionMatrix).replaceAll("[\\[\\],]", ""));
		System.out.println(Arrays.deepToString(emmissionMatrix).replaceAll("[\\[\\],]", ""));
		}
		System.out.println(prob_model);
		
	}
	
	public static void main(String[] args){
		String[] input =null;
		BaumWelch calc = new BaumWelch();
		
		try {
			input =  calc.readInput();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
//		for(int i=0;i<4;i++){
//			System.out.println(input[i]);
//		}
		
		calc.calculateEmissionProbability(input);
/*		System.out.println(Arrays.deepToString(calc.alpha));
		System.out.println(Arrays.deepToString(calc.beta));		
		System.out.println("===");
		System.out.println(Arrays.deepToString(calc.digamma));
		System.out.println(Arrays.deepToString(calc.gamma));
		System.out.println(calc.T);
		System.out.println(calc.N);
		System.out.println(calc.K);*/
		//System.out.println(calc.prob_model);

		
	}
	

}
