import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Scanner;

/**
 * This class calculates the next emission distribution
 */

/**
 * @author akash
 *
 */
public class NextEmissionPredictor {
	
	private double[][] transitionMatrix = null;
	private double[][] emmissionMatrix = null;
	private double[][] stateDistribution = null;
	private double[][] nextEmissionDistribution = null;
	private double[][] nextStateDistribution = null;
	/*
	 * reads input file through pipe
	 */
	private  String[] readInput() throws IOException{
		String[] input = new String[3];
		InputStreamReader inputStream = new InputStreamReader(System.in);
		BufferedReader bufferedReader = new BufferedReader(inputStream);
		int i =0;
		while(i<3){
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
		
	private  double[][] calculateNextEmissionDistribution(String[] hmmModel){
		transitionMatrix = convertStringTo2DArray(hmmModel[0]);
		emmissionMatrix = convertStringTo2DArray(hmmModel[1]);
		stateDistribution = convertStringTo2DArray(hmmModel[2]);
		nextStateDistribution = matrixMultiplication(stateDistribution, transitionMatrix);
		nextEmissionDistribution = matrixMultiplication(nextStateDistribution,emmissionMatrix);
		return nextEmissionDistribution;
	}
	
	public static void main(String[] args){
		String[] input =null;
		NextEmissionPredictor predictor1 = new NextEmissionPredictor();
		
		try {
			input =  predictor1.readInput();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//uncomment for debugging
		/*for(int i=0;i<3;i++){
			System.out.println(input[i]);
		}*/
		
		double[][] nextEmissionProbabilities = predictor1.calculateNextEmissionDistribution(input);
		String output = "1"+" "+nextEmissionProbabilities[0].length+" "+Arrays.deepToString(nextEmissionProbabilities).replaceAll("[\\[\\],]", "");
		System.out.println(output);
	}

	
}