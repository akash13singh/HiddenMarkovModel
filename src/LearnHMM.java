import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

/**
 * Code to calculate emission sequence probability. Implements the Forward Pass Algorithm
 * Refer Stamp Tutorial: https://www.cs.sjsu.edu/~stamp/RUA/HMM.pdf
 */

/**
 * @author akash
 *
 */
public class LearnHMM {

	/*
	 * TODO use setters and getters for these variables.
	 */
	private double[][] transitionMatrix = null;
	private double[][] emissionMatrix = null;
	private double[][] initialStateDistribution = null;
	private int[] emissionSequence = null;
	
	// number of states
	private int N;
	
	// length of emmission sequence;
	private int T;
	
	//no. of possible emissions
	private int M;
	
	//c: scale factor
	private double[] c = null;
	
	// an N*T matrix used in alpha/forward pass
	private double[][] alpha = null;
	
	private double[][] beta = null;
	
	private double[][] gamma = null;
	
	private double[][][] digamma = null;

	
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
	
	
	private void alphaPass(){
		
	    
		c[0]=0;
	    //calculate alpha 0;
		for(int i = 0; i<N; i++){
			alpha[i][0]= initialStateDistribution[0][i] * emissionMatrix[i][emissionSequence[0]];
			c[0] = c[0] + alpha[i][0];
		}		
		
		//scale alpha[i][0]
		c[0] = 1.0/c[0];	
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
				alpha[i][t] = alpha[i][t] * emissionMatrix[i][emissionSequence[t]];
				c[t]= c[t]+alpha[i][t];
			}
			
			//scale alpha[i][t]
			c[t] = 1.0/c[t];
			for(int i = 0; i<N; i++){
				alpha[i][t] = c[t]*alpha[i][t];	
			}
		}	
	}
	
	private void betaPass(){
		
		//initialize beta[T-1]
		for(int i = 0; i<N; i++){
			beta[i][T-1]= c[T-1];
		}	
		
		for(int t = T-2;t>=0;t--){
			for(int i =0; i<N;i++){
				beta[i][t]=0;
				for(int j =0;j<N;j++){
					beta[i][t] = beta[i][t] + (transitionMatrix[i][j] * 
							     emissionMatrix[j][emissionSequence[t+1]] * beta[j][t+1]);
				}
				beta[i][t] = c[t] * beta[i][t];
			}		
		}		
	}
	
	private void computeGammas(){
		
		for(int t=0;t<=T-2;t++){
			
			double denom = 0.0;
			for(int i =0;i<N;i++){
				for(int j=0;j<N;j++){
					denom = denom + (alpha[i][t] * transitionMatrix[i][j]* emissionMatrix[j][emissionSequence[t+1]] * beta[j][t+1]);
				}
			}
			
			for(int i =0;i<N;i++){
				gamma[i][t] = 0.0;
				for(int j=0;j<N;j++){
					digamma[i][j][t] = (alpha[i][t] * transitionMatrix[i][j]* emissionMatrix[j][emissionSequence[t+1]] * beta[j][t+1])/denom;
					gamma[i][t] = gamma[i][t] + digamma[i][j][t];
				}
			}			
		}
		
		//for gamma[T-1]
		double denom = 0;
		for(int i =0;i<N;i++){
			denom = denom + alpha[i][T-1];
		}
		for(int i =0;i<=N-1;i++){
			gamma[i][T-1] = alpha[i][T-1]/denom;
		}
	}
	
	private void reEstimate(){
		
		//reEstimate pi
		for(int i =0;i<N;i++){
			initialStateDistribution[0][i] = gamma[i][0];
		}
		
		//re-estimate A
		for(int i =0;i<N;i++){
			for(int j =0;j<N;j++){
				double numer =0.0;
				double denom=0.0;
				for(int t =0;t<=T-2;t++){
					numer = numer + digamma[i][j][t];							
					denom = denom + gamma[i][t];
				}
				transitionMatrix[i][j] = numer/denom;
			}
		}
		
		//re-estimate B
		for(int i =0;i<N;i++){
			for(int j =0;j<M;j++){
				double numer=0.0;
				double denom=0.0;
				for(int t =0;t<T;t++){
					if(emissionSequence[t]==j) {
						numer = numer + gamma[i][t];
					}
					denom = denom + gamma[i][t];
				}
				emissionMatrix[i][j] = numer/denom;
			}
		}
		
	}
	
	
	private void braumWelch(){
		
		int maxIters =50;
		int iters = 0;
		double oldLogProb = Double.NEGATIVE_INFINITY;
		
		while(true){
			alphaPass();
			betaPass();
			computeGammas();
			reEstimate();
			double logProb = 0;
			
			for(int i =0;i<T;i++){
				logProb = logProb + Math.log(c[i]);
			}
			
			logProb = -logProb;
			iters=iters+1;
			if(iters<maxIters /*&& logProb > oldLogProb*/){
				oldLogProb = logProb;
			}else{
				break;
			}
		}
		
		//System.out.println("Iterations = "+iters);
		for(int i=0;i<N;i++){
			for(int j=0;j<N-1;j++){
				transitionMatrix[i][j] = new BigDecimal(transitionMatrix[i][j]).setScale(6, RoundingMode.HALF_UP).doubleValue();
			}
			for(int m=0;m<M-1;m++){
				emissionMatrix[i][m] = new BigDecimal(emissionMatrix[i][m]).setScale(6, RoundingMode.HALF_UP).doubleValue();
			}
		}
		
		System.out.println(N+" "+N+" "+Arrays.deepToString(transitionMatrix).replaceAll("[\\[\\],]", ""));
		System.out.println(N+" "+M+" "+Arrays.deepToString(emissionMatrix).replaceAll("[\\[\\],]", ""));
	}
	
	
	private void initialize(String[] hmmParams){	
		transitionMatrix = convertStringTo2DArray(hmmParams[0]);
		emissionMatrix = convertStringTo2DArray(hmmParams[1]);
		initialStateDistribution = convertStringTo2DArray(hmmParams[2]);
		emissionSequence = convertStringTo1DArray(hmmParams[3]);
		
		
		N = transitionMatrix.length;
		T= emissionSequence.length;	
		M = emissionMatrix[0].length;
		
		c = new double[T];
		alpha = new double[N][T];
		beta = new double[N][T];
		gamma = new double[N][T];
		digamma = new double[N][N][T];

	}
	
	public static void main(String[] args){
		String[] input =null;
		LearnHMM hmm = new LearnHMM();
		
		try {
			input =  hmm.readInput();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
//		for(int i=0;i<4;i++){
//			System.out.println(input[i]);
//		}
		
		hmm.initialize(input);
		hmm.braumWelch();
		
		
	}
	

}
