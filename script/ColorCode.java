import java.util.*;
import java.io.*;
public class ColorCode{
  public static void main(String argv[]){
    new ColorCode(argv);
  }
  public ColorCode(String argv[]){
    if(argv.length<2){
      System.out.println("polyfile objfile");
      System.exit(0);
    }
  //off is zero based. obj is 1 based
    Scanner polysc = openFile(argv[0]);
    Scanner objsc = openFile(argv[1]); 
    //build polygons as sets of vertex indices
    //maps from a vertex to all its adjacent polygons
    ArrayList<Integer> v2poly[] ;
    HashSet<Integer>[] poly;
    polysc.nextLine();
    int nVert, nPoly, nTrig=0;
    String line = polysc.nextLine();
    String toks[] = line.split(" ");
    nVert = Integer.parseInt(toks[0]);
    nPoly = Integer.parseInt(toks[1]);
    v2poly = (ArrayList<Integer>[])new ArrayList[nVert];
    poly = (HashSet<Integer>[])new HashSet[nPoly];
    polysc.nextLine();
    for(int ii=0;ii<nVert;ii++){
      polysc.nextLine();
    }
    for(int ii=0;ii<nPoly;ii++){
        line = polysc.nextLine();
        toks = line.split("\\s+");
        HashSet<Integer> s  =new HashSet<Integer>();
        for(int jj=1;jj<toks.length;jj++){
            int idx = Integer.parseInt(toks[jj]);
            s.add(idx);
            if(v2poly[idx]==null){
                v2poly[idx] = new ArrayList<Integer>();
            }
            v2poly[idx].add(ii);
        }
        nTrig += toks.length-3;
        poly[ii] = s;
    }
    polysc.close();
    System.out.println(nVert);
    System.out.println(nTrig);
    while(objsc.hasNextLine()){
        line = objsc.nextLine();
        if(line.length()<1){
            continue;
        }
        if(line.charAt(0)=='#'){
            continue;
        }
        if(line.charAt(0)=='f'){
            toks = line.split(" ");
            int idx[]= new int[3];
            HashSet<Integer> verts = new HashSet<Integer>();
            for(int ii=1;ii<toks.length;ii++){
                idx[ii-1] = Integer.parseInt(toks[ii])-1;
                verts.add(idx[ii-1]);
            }
            ArrayList<Integer> polylist = v2poly[idx[0]];
            int polyIdx=0;
            for(int ii=0;ii<polylist.size();ii++){
                int kk = polylist.get(ii);
                if(poly[kk].containsAll(verts)){
                    polyIdx = kk;
                    break;
                }
            }
            System.out.println("3 "+idx[0]+" "+idx[1]+" "+idx[2]+" "+polyIdx);
        }
        else if(line.charAt(0)=='v'){
            toks = line.split(" ");
            System.out.println(toks[1]+" "+toks[2]+" "+toks[3]);
        }
            
    }
    objsc.close();
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
