import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
<<<<<<< HEAD
 * Code to calculate emission sequence probability. Implements the Forward Pass Algorithm
 * Refer Stamp Tutorial: https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
 */

/**
 * @author mazen
=======
 * Implements the Viterbi pass algorith. A DP algorithm to calculate the most likely state sequence 
 * given the HMM Model and observation sequence (emission sequence).
 */

/**
 * @author akash
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
 *
 */
public class Viterbi {

	/*
	 * TODO use setters and getters for these variables.
	 */
	private double[][] transitionMatrix = null;
	private double[][] emmissionMatrix = null;
	private double[][] initialStateDistribution = null;
	private int[] emissionSequence = null;
<<<<<<< HEAD
	//c: scale factor
	private double[] c = null;
	
	private int[] stateSeq = null;
	private double[] stateSeqProb = null;

	
	
	// an N*T matrix used in alpha/forward pass
	private double[][] alpha = null;
	
=======
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
	
	// number of states
	private int N;
	
	// length of emmission sequence;
	private int T;
	
<<<<<<< HEAD
=======
	//N*T matrix to store δt(i) the prob of having observed O[1:t] and being in state i given the most likely state for each t
	private double[][] delta = null;
	
	//N*T matrix to store the index of the most likely state
	private int[][] deltaIdx = null;
	
	// a T length array to store the most likely state sequnce
	private int[] mostLikelyStates = null;

	
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
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
<<<<<<< HEAD
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
	
	private void viterbi(){
		
	    c = new double[T];
		alpha = new double[N][T];
		stateSeq = new int[T];
		stateSeqProb = new double[T];

	    //calculate alpha 0;
		//c[0]=0 ;
		for(int i = 0; i<N; i++){
			alpha[i][0]= initialStateDistribution[0][i] * emmissionMatrix[i][emissionSequence[0]];
			c[0] = c[0] + alpha[i][0];
		}		
		c[0] = 1/c[0];	
		//System.out.println(Arrays.deepToString(alpha));

		//scale alpha[i][0]
		for(int i = 0; i<N; i++){
			alpha[i][0]= c[0] * alpha[i][0]; 
		}
		
		double max = -999;
		int max_state = -1 ;
		for(int i = 0; i<N; i++){
			if ( alpha[i][0] > max  ){
				max = alpha[i][0] ; 
				max_state = i;
				}
		}
		stateSeq[0] = max_state;
		stateSeqProb[0] = max;
		//System.out.println(max);
		
		//System.out.println(Arrays.toString(stateSeq));

		
		//compute alpha[i][t]
		for(int t =1; t<T; t++){
			c[t]=0;
			for(int i = 0; i<N; i++){
				//alpha[i][t]= stateSeqProb[t-1] * transitionMatrix[stateSeq[t-1] ][i];
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
/*			System.out.println(Arrays.deepToString(alpha));
			System.out.println("============");*/

			max = -999;
		    max_state = -1 ;
			for(int i = 0; i<N; i++){
				if ( alpha[i][t] > max  ){
					max = alpha[i][t] ; 
					max_state = i;
					}
			}
			stateSeq[t] = max_state;
			stateSeqProb[t] = max;
			//System.out.println(Arrays.toString(stateSeq));
		}
		
		
	}
	
	private double calculateEmissionProbability(String[] hmmParams){
		double prob =1;
=======
	 * implememts the viterbi DP logic.
	 * uses log to store delta values to avoid underflow.
	 */
	private double viterbiDP(){
		
		delta = new double[N][T];
		deltaIdx = new int[N][T];
		mostLikelyStates = new int[T];
		
	    //calculate delta 0;
		for(int i = 0; i<N; i++){
			//delta[i][0]= Math.log(initialStateDistribution[0][i]) + Math.log( emmissionMatrix[i][emissionSequence[0]]);
			delta[i][0]= (initialStateDistribution[0][i] * emmissionMatrix[i][emissionSequence[0]]);
		}
		
		//compute delta[i][t]
		for(int t =1; t<T; t++){
			for(int i = 0; i<N; i++){
				delta[i][t]=0;
				
				double max = -1000;
				int idx=-1;
				//double logEmissinProb = Math.log(emmissionMatrix[i][emissionSequence[t]]);
				
				for(int j=0; j<N;j++){

					double tmp =  (delta[j][t-1] *transitionMatrix[j][i] * emmissionMatrix[i][emissionSequence[t]]) ;
					if(tmp>max){
						max=tmp;
						idx = j;
					}
				}
				delta[i][t] = max;
				deltaIdx[i][t] = idx;
				mostLikelyStates[t] =idx;
			}
		}
		
		//System.out.println(Arrays.deepToString(delta));
		
		//find prob of max likely sequence;
		double maxProb =0;
		int maxIdx = 0;
		
//		for(int i =0 ;i<N;i++){
//			for(int j =0 ;j<T;j++){
//				System.out.print(delta[i][j]+" , ");
//			}
//			System.out.println();
//		}
//
//		for(int i =0 ;i<N;i++){
//			for(int j =0 ;j<T;j++){
//				System.out.print(deltaIdx[i][j]+" , ");
//			}
//			System.out.println();
//		}
		
		for(int i=0; i<N;i++){
			if(delta[i][T-1] > maxProb){
				maxProb  = delta[i][T-1];
				maxIdx = i;
			}
		}
		
		mostLikelyStates[T-1] = maxIdx;	
		for(int t =T-2;t>=0;t--){
			mostLikelyStates[t]= deltaIdx[mostLikelyStates[t+1]][t+1];
		}
			
		return maxProb;
		
	}
	
	private void estimateLikelyStateSEquence(String[] hmmParams){

>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
		transitionMatrix = convertStringTo2DArray(hmmParams[0]);
		emmissionMatrix = convertStringTo2DArray(hmmParams[1]);
		initialStateDistribution = convertStringTo2DArray(hmmParams[2]);
		emissionSequence = convertStringTo1DArray(hmmParams[3]);
		
<<<<<<< HEAD
//		uncomment for debugging
=======
		//uncomment for debugging
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
//		System.out.println(Arrays.deepToString(transitionMatrix));
//		System.out.println(Arrays.deepToString(emmissionMatrix));
//		System.out.println(Arrays.deepToString(initialStateDistribution));
//		System.out.println(Arrays.toString(emissionSequence));
		
		N = transitionMatrix.length;
		T= emissionSequence.length;
		
<<<<<<< HEAD
		viterbi();
		
		for(int t = 0; t<T; t++){
			prob = prob*c[t];
		}
		return 1/prob;
=======
		 double maxProb = viterbiDP();
		 //System.out.println("Max Prob :"+maxProb);
		 System.out.println(Arrays.toString(mostLikelyStates).replaceAll("[\\[\\],]", ""));
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
	}
	
	public static void main(String[] args){
		String[] input =null;
<<<<<<< HEAD
		Viterbi vit = new Viterbi();
		
		try {
			input =  vit.readInput();
=======
		Viterbi vitrb = new Viterbi();
		
		try {
			input =  vitrb.readInput();
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
//		for(int i=0;i<4;i++){
//			System.out.println(input[i]);
//		}
		
<<<<<<< HEAD
		double prob = vit.calculateEmissionProbability(input);
		System.out.println(Arrays.toString(vit.stateSeq).replaceAll("[\\[\\],]", ""));
=======
		vitrb.estimateLikelyStateSEquence(input);
		
>>>>>>> 366b6870d045c02917ac3f75602397f4554435eb
	}
	

}
