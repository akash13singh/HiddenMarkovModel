import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

/**
 * Code to calculate emission sequence probability. Implements the Forward Pass Algorithm
 * Refer Stamp Tutorial: https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
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
	//c: scale factor
	private double[] c = null;
	
	private int[] stateSeq = null;
	
	// an N*T matrix used in alpha/forward pass
	private double[][] alpha = null;
	
	
	// number of states
	private int N;
	
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
		stateSeq = new int[T];
		
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
		
		double max = -999;
		int max_state = -1 ;
		for(int i = 0; i<N; i++){
			if ( alpha[i][0] > max  ){
				max = alpha[i][0] ; 
				max_state = i;
				}
		}
		
/*		System.out.println(max);
		stateSeq[0] = max_state;
		System.out.println(Arrays.toString(stateSeq));*/

		
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
			System.out.println(Arrays.deepToString(alpha));
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
			//System.out.println(Arrays.toString(stateSeq));
		}
		
		
	}
	
	private double calculateEmissionProbability(String[] hmmParams){
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
		T= emissionSequence.length;
		
		alphaPass();
		
		for(int t = 0; t<T; t++){
			prob = prob*c[t];
		}
		return 1/prob;
	}
	
	public static void main(String[] args){
		String[] input =null;
		Viterbi vit = new Viterbi();
		
		try {
			input =  vit.readInput();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
//		for(int i=0;i<4;i++){
//			System.out.println(input[i]);
//		}
		
		double prob = vit.calculateEmissionProbability(input);
		System.out.println(Arrays.toString(vit.stateSeq).replaceAll("[\\[\\],]", ""));
		
	}
	

}
