import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
 * Implements the Viterbi pass algorith. A DP algorithm to calculate the most likely state sequence 
 * given the HMM Model and observation sequence (emission sequence).
 */

/**
 * @author akash
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
	
	// number of states
	private int N;
	
	// length of emmission sequence;
	private int T;
	
	//N*T matrix to store δt(i) the prob of having observed O[1:t] and being in state i given the most likely state for each t
	private double[][] delta = null;
	
	//N*T matrix to store the index of the most likely state
	private int[][] deltaIdx = null;
	
	// a T length array to store the most likely state sequnce
	private int[] mostLikelyStates = null;

	
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

		transitionMatrix = convertStringTo2DArray(hmmParams[0]);
		emmissionMatrix = convertStringTo2DArray(hmmParams[1]);
		initialStateDistribution = convertStringTo2DArray(hmmParams[2]);
		emissionSequence = convertStringTo1DArray(hmmParams[3]);
		
		//uncomment for debugging
//		System.out.println(Arrays.deepToString(transitionMatrix));
//		System.out.println(Arrays.deepToString(emmissionMatrix));
//		System.out.println(Arrays.deepToString(initialStateDistribution));
//		System.out.println(Arrays.toString(emissionSequence));
		
		N = transitionMatrix.length;
		T= emissionSequence.length;
		
		 double maxProb = viterbiDP();
		 //System.out.println("Max Prob :"+maxProb);
		 System.out.println(Arrays.toString(mostLikelyStates).replaceAll("[\\[\\],]", ""));
	}
	
	public static void main(String[] args){
		String[] input =null;
		Viterbi vitrb = new Viterbi();
		
		try {
			input =  vitrb.readInput();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
//		for(int i=0;i<4;i++){
//			System.out.println(input[i]);
//		}
		
		vitrb.estimateLikelyStateSEquence(input);
		
	}
	

}