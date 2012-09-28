import java.util.*;
import java.io.*;
public class ColorQSlim{
  public static void main(String argv[]){
    new ColorQSlim(argv);
  }
  public ColorQSlim(String argv[]){
    if(argv.length<1){
      System.out.println("objfile");
      System.exit(0);
    }
  // obj is 1 based
    Scanner objsc = openFile(argv[1]);
    ArrayList<String[]> faces;
    ArrayList<String[]> verts;
    while(objsc.hasNextLine()){
        line = objsc.nextLine();
        if(line.length()<1){
            continue;
        }
        if(line.charAt(0)=='#'){
            continue;
        }
        if(line.charAt(0)=='f'){
            faces.add(line);
        }
        else if(line.charAt(0)=='v'){
            verts.add(line);
        }
            
    }
    objsc.close();
  
    System.out.println(verts.size(););
    System.out.println(faces.size());
   
    }
  Scanner openFile(String filename){
    Scanner sc = null;
    try{
      sc = new Scanner (new FileInputStream(filename));
    }catch(IOException ioe){
      System.out.println(ioe.getMessage());
    }
    return sc;
  }
}
