import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.io.*; 
import javax.crypto.spec.PSource;
// this program creates windowed feature files for single cysteine
public class run_windowing{
	
	public static void main(String[] args) throws IOException {	
		
		
		
		// BufferedReader brTest = new BufferedReader(new FileReader("../../../Input/Dataset_Name.txt"));
 	// 	String Dataset_Name=brTest .readLine();
 	// 	System.out.println(Dataset_Name);

		BufferedReader inputs = new BufferedReader(new FileReader("./window_param.txt"));
		int num_of_features=Integer.parseInt(inputs.readLine());
		int start_ws=Integer.parseInt(inputs.readLine());
		int end_ws=Integer.parseInt(inputs.readLine());
		// System.out.println(Dataset_Name);
 		System.out.println(num_of_features);
 		System.out.println(start_ws);
 		System.out.println(end_ws);

 		// int num_of_features = 112; // you should provide exact number of features here without Target
	
		String idListPath = "../Features/id_list.txt";
		String combinedFeatureDir = "../output/merge_features/";
		String singleCysWindowedFilesDir = "../output/Windowed_file/";
		
		for(int ws = start_ws; ws <= end_ws; ws=ws+2){
			
			String singleCysWindowedFileName = singleCysWindowedFilesDir+"feat_"+num_of_features+"_w_"+ws+".csv";
			System.out.println(singleCysWindowedFileName);
			createWindow(idListPath, combinedFeatureDir, singleCysWindowedFileName, ws, num_of_features);			
			
		}
		
	}
	public static void createWindow(String idListPath, String combinedFeatureDir, String singleCysWindowedFileName, int window_size, int number_of_features) throws IOException {
			
		int temp=0;
		BufferedWriter windowedFileWriter = BufferReaderAndWriter
				.getWriter(new File(singleCysWindowedFileName));
				//??
		
		// windowedFileWriter.write("class,");
		
		// for(int i = 0; i < 525; i++){
			
		// 	if(i == 524){
				
		// 		windowedFileWriter.write("V"+(i+1));
				
		// 	}else{
				
		// 		windowedFileWriter.write("V"+(i+1)+",");
				
		// 	}			
			
		// }
		
		// windowedFileWriter.write("\n");
				//??
		
		// parse the list file which contains protein ids
		
		BufferedReader listFileRd = BufferReaderAndWriter.getReader(new File(idListPath));
		String proteinID = "";
		
	
		while ((proteinID = listFileRd.readLine())!=null) {
	
			// loading file in memory
			List<String> featuresArray = new ArrayList<String>();
			
			System.out.println(proteinID);
				
			String featureFilePath = combinedFeatureDir+proteinID+".csv";
			
			BufferedReader featureFileReader = BufferReaderAndWriter.getReader(new File(featureFilePath));
			
			String featureLine = "";
			
			while((featureLine = featureFileReader.readLine())!=null){
			
				String[] element = featureLine.split(",");
				if(element[0].equalsIgnoreCase("Target")){
					
					continue;					
					
				}
				
				featuresArray.add(featureLine);			
				temp++;
			}
			
			featureFileReader.close();
			
						
			// start creating features based on window size for each residue
			
			int beforeAndAfterRowsToBeIncluded = window_size/2;
			
			
			for (int j = 0; j < featuresArray.size(); j++) {
				
				String featuresString = "";
			//	int featureCounter = 0;
				
				
	//			/*=================== Ignore the cystein pairs, either of which can not form consecutive residues equal to the window size =============*/
	//			int beforeWindowCheck = j - beforeAndAfterRowsToBeIncluded;
	//			int afterWindowCheck = j + beforeAndAfterRowsToBeIncluded;
	//			
	//			if((beforeWindowCheck < 0) || (afterWindowCheck >= featuresArray.size())){
	//				System.out.println(proteinID+"\t"+featuresArray.size()+"\t"+j+"\t"+j);
	//				continue;					
	//				
	//			}
				
				if(beforeAndAfterRowsToBeIncluded > 0){		// if window size is > 1
					
					int beforeFeatureCounter = beforeAndAfterRowsToBeIncluded;
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the left
																		
						int beforeIndex = j-beforeFeatureCounter;
						
						if(beforeIndex < 0){										// if for terminal residues features does not exist
							
							for(int l = 0; l < number_of_features; l++){
							//	featureCounter++;
//								featuresString = featuresString+"9999.9999,";
								featuresString = featuresString+"0,";
								
							}
							
						}else {
							
							String[] beforeFeatureExist = featuresArray.get(beforeIndex).trim().split(",");
							
							for(int m = 1; m < beforeFeatureExist.length; m++){
							//	featureCounter++;
							// 	featuresString = featuresString+featureCounter+":"+beforeFeatureExist[m]+" ";
							//	featuresString = featuresString+beforeFeatureExist[m]+" ";
	//							String[] feature = beforeFeatureExist[m].split(":");
								
								featuresString = featuresString+beforeFeatureExist[m]+",";
								
								
							}
							
							
						}
						
						beforeFeatureCounter--;
												
					}
					
					String[] currentIndexFeature = featuresArray.get(j).trim().split(",");
					
					for(int n = 1; n < currentIndexFeature.length; n++){
					//	featureCounter++;
					//	featuresString = featuresString+featureCounter+":"+currentIndexFeature[n]+" ";
					//	featuresString = featuresString+currentIndexFeature[n]+" ";
						
	//					String[] feature = currentIndexFeature[n].split(":");
						featuresString = featuresString+currentIndexFeature[n]+",";
						
					}
					
	//				featuresString = currentIndexFeature[0]+","+featuresString;
					
										
					int afterFeatureCounter = 0;
					
					for(int k = 0; k < beforeAndAfterRowsToBeIncluded; k++){		// append the features to the right
						
						afterFeatureCounter++;
						
						int afterIndex = j+afterFeatureCounter;
						
						if(afterIndex > featuresArray.size()-1){
							
							for(int l = 0; l < number_of_features; l++){
							//	featureCounter++;
							//	featuresString = featuresString+featureCounter+":0 ";
								
//								featuresString = featuresString+"9999.9999,";
								featuresString = featuresString+"0,";
								
							}
							
						}else {
							
							String[] afterFeatureExist = featuresArray.get(afterIndex).trim().split(",");
							
							for(int m = 1; m < afterFeatureExist.length; m++){
							//	featureCounter++;
							//	featuresString = featuresString+featureCounter+":"+afterFeatureExist[m]+" ";
							//	featuresString = featuresString+afterFeatureExist[m]+" ";
	//							String[] feature = afterFeatureExist[m].split(":");
								
								featuresString = featuresString+afterFeatureExist[m]+",";
								
							}
							
																					
						}
												
					}
					
					featuresString = currentIndexFeature[0]+","+featuresString;
					
					
				}else { // if window size is 1
					
					String[] thisFeature = featuresArray.get(j).trim().split(",");
					
					for(int p = 0; p < thisFeature.length; p++){
						
					//	featureCounter++;
					//	featuresString = featuresString+featureCounter+":"+thisFeature[p]+" ";
					//	featuresString = featuresString+thisFeature[p]+" ";
	//					String[] feature = thisFeature[p].split(":");
						
						featuresString = featuresString+thisFeature[p]+",";							
						
					}
					
	//				featuresString = thisFeature[0]+","+featuresString.trim();
					
										
				}
				
				if(featuresString.contains("inf") || featuresString.contains("nan")){
					
					System.out.println("protein contains inf "+proteinID);
					System.exit(0);
				}
				//System.out.println(featuresString.substring(0, featuresString.length()-1));
				//windowedFileWriter.write(proteinID);
				windowedFileWriter.write(featuresString.substring(0, featuresString.length()-1));
//				windowedFileWriter.write(featuresString);
				windowedFileWriter.write("\n");
				windowedFileWriter.flush();
				
			}
			
		}
		System.out.println(temp);
		listFileRd.close();
		windowedFileWriter.close();	
	
	}
}
