import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * 
 */

/**
 * @author akash
 *
 */
public class NextEmissionProbablityCalculator {

	private double[][] transitionMatrix = null;
	private double[][] emmissionMatrix = null;
	private double[][] initialStateDistribution = null;
	private double[][] emissionSequence = null;
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
	

}
