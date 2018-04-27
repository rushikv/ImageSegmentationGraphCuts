/******************************************************************************
 *  Compilation:  javac UI.java
 *  Execution:    java UI image_file_name
 *  Dependencies: CountMinCuts.java
 *
 *  UI for image segmentation using Rachel Silva's algorithm
 *
 * Author: Rushik Vartak (rv9981@g.rit.edu)
 *
 ******************************************************************************/
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.*;
import java.util.ArrayList;

import javax.imageio.ImageIO;
import javax.swing.*;



public class UI extends JPanel implements MouseListener, MouseMotionListener{
    int mouse_x, mouse_y, x,y;
    static int x_prev = -1;
    static int y_prev = -1;
    static BufferedImage image;
    static BufferedImage imageOG;
    Dimension size = new Dimension();
    String imgFile;
    static int count=0;
    static adj_list G;
    static ArrayList<Integer> S = new ArrayList<>(){};
    static ArrayList<Integer> T = new ArrayList<>(){};
    boolean scrib_s, scrib_t;
    Menu menu;
    Menu menu2;
    Menu menu3;
    static ArrayList<int[][]> finalCuts;
    static int type, samples;
    static int pixCounts[][];



    public UI(BufferedImage image) {
        this.image = image;
        size.setSize(image.getWidth(), image.getHeight());
        this.addMouseListener(this);
        this.addMouseMotionListener(this);
    }


    protected void paintComponent(Graphics g){
        super.paintComponent(g);
        int x = (getWidth() - size.width)/2;
        int y = (getHeight() - size.height)/2;
        g.drawImage(image, x, y, this);
        g.setColor(Color.blue);
        g.drawLine(mouse_x,mouse_y-15,mouse_x,mouse_y-5);
        g.drawLine(mouse_x-5,mouse_y-10,mouse_x+5,mouse_y-10);
    }

    public static void main(String[] args){

        SwingUtilities.invokeLater(new Runnable()
        {
            public void run()
            {
                new UI();

            }
        });
    }

    public UI() {
        try{
            initialize();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    private void initialize(){
        loadImage();
    }

    public void loadImage() {
//        imgFile = "295087.jpg";
//        imgFile = "imgs/118.jpg";
//        imgFile = "test_img3.png";
//        imgFile = "testsame2.png";
//        imgFile = "oscartest.png";
//        imgFile = "test10x10.png";
        imgFile = "perf_1200x1200.jpg";
        try {
            image = ImageIO.read(new File(imgFile));
            imageOG = ImageIO.read(new File(imgFile));
        } catch (IOException e) {
            e.printStackTrace();
        }
        UI test = new UI(image);
        JFrame f = new JFrame();

        f.setTitle("Viewer");
        f.add(new JScrollPane(test));
        f.setVisible(true);
        Insets insets = f.getInsets();
        f.setSize(image.getWidth()+insets.left+insets.right+1, image.getHeight()+insets.top+insets.bottom+1);

        f.setResizable(false);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        MenuBar menuBar = new MenuBar();
        menu = new Menu("Scribble");
        menu.add("Draw_S");
        menu.add("Draw_T");
        menu.getItem(1).setEnabled(false);

        menu3 = new Menu("Compute");
        menu3.add("Left_Right_Most");
        menu3.add("Prob_Sample");
        menu3.setEnabled(false);

        menu2 = new Menu("Filter");
        menu2.add("Grayscale");
        menu2.add("RGB");
        menu2.setEnabled(false);

        count = 0;
        menuBar.add(menu);
        menuBar.add(menu3);
        menuBar.add(menu2);

        menu.addActionListener(new axnListener());
        menu2.addActionListener(new axnListener2());
        menu3.addActionListener(new axnListener3());
        f.setMenuBar(menuBar);
    }

    public class axnListener implements ActionListener{
        public void actionPerformed(ActionEvent e){
            if(e.getActionCommand().equalsIgnoreCase("Draw_S")){
                scrib_s = true;
                scrib_t = false;
                count++;
                menu.getItem(0).setEnabled(false);
                menu.getItem(1).setEnabled(true);
                menu3.setEnabled(false);
                menu2.setEnabled(false);

            }
            if(e.getActionCommand().equalsIgnoreCase("Draw_T")){
                scrib_t = true;
                scrib_s = false;
                x_prev = -1;
                y_prev = -1;
                count++;
                menu.setEnabled(false);
                menu3.setEnabled(true);
                menu2.setEnabled(false);

            }
        }
    }

    public class axnListener2 implements ActionListener{
        public void actionPerformed(ActionEvent e){
            if(e.getActionCommand().equalsIgnoreCase("Grayscale")){

                for(int xi=0;xi<image.getWidth();xi++){
                    for(int yi=0;yi<image.getHeight();yi++){
                        image.setRGB(xi,yi,RGBInt((int)luminence(RGB(image.getRGB(xi,yi)))));
                    }
                }

                if(type==0) {
                    int[][] finalCut = finalCuts.get(0);
                    for (int i = 0; i < finalCut.length; i++) {
                        int[] p1 = vertexIJ(finalCut[i][0]);
                        int[] p2 = vertexIJ(finalCut[i][1]);
                        image.setRGB(p1[0], p1[1], 255);
                        image.setRGB(p2[0], p2[1], 255);
                    }
                    if (finalCuts.size() > 1) {
                        finalCut = finalCuts.get(1);
                        for (int i = 0; i < finalCut.length; i++) {
                            int[] p1 = vertexIJ(finalCut[i][0]);
                            int[] p2 = vertexIJ(finalCut[i][1]);
                            image.setRGB(p1[0], p1[1], 16776960);
                            image.setRGB(p2[0], p2[1], 16776960);
                        }
                    }
                }
                else if(type==1){
                    for(int xi=0;xi<image.getWidth();xi++){
                        for(int yi=0;yi<image.getHeight();yi++) {
                            if(pixCounts[xi][yi] != 0){
                                int value = (int)(((double)pixCounts[xi][yi]/samples)*256);
                                image.setRGB(xi,yi,16776960+(256-value));
                            }
                        }
                    }
                }
                System.out.println("------------------x-------------------------x-----------------------------");

            }
            else if(e.getActionCommand().equalsIgnoreCase("RGB")){

                for(int xi=0;xi<image.getWidth();xi++){
                    for(int yi=0;yi<image.getHeight();yi++){
                        image.setRGB(xi,yi,imageOG.getRGB(xi,yi));
                    }
                }

                if(type==0) {
                    int[][] finalCut = finalCuts.get(0);
                    for (int i = 0; i < finalCut.length; i++) {
                        int[] p1 = vertexIJ(finalCut[i][0]);
                        int[] p2 = vertexIJ(finalCut[i][1]);
                        image.setRGB(p1[0], p1[1], 255);
                        image.setRGB(p2[0], p2[1], 255);
                    }
                    if (finalCuts.size() > 1) {
                        finalCut = finalCuts.get(1);
                        for (int i = 0; i < finalCut.length; i++) {
                            int[] p1 = vertexIJ(finalCut[i][0]);
                            int[] p2 = vertexIJ(finalCut[i][1]);
                            image.setRGB(p1[0], p1[1], 16776960);
                            image.setRGB(p2[0], p2[1], 16776960);
                        }
                    }
                }
                else if(type==1){
                    for(int xi=0;xi<image.getWidth();xi++){
                        for(int yi=0;yi<image.getHeight();yi++) {
                            if(pixCounts[xi][yi] != 0){
                                int value = (int)(((double)pixCounts[xi][yi]/samples)*256);
                                image.setRGB(xi,yi,16776960+(256-value));
                            }
                        }
                    }
                }
                System.out.println("------------------x-------------------------x-----------------------------");

            }
        }
    }

    public class axnListener3 implements ActionListener{
        public void actionPerformed(ActionEvent e){
            if(e.getActionCommand().equalsIgnoreCase("Left_Right_Most")){

                type = 0;
                menu2.setEnabled(true);
                menu.setEnabled(false);


                count++;
                scrib_s = false;
                scrib_t = false;

//                S.add(vertexNo(25,19));
//                T.add(vertexNo(17,19));

                System.out.print("Source S: [");
                for(int ss=0;ss<S.size();ss++){
                    System.out.print("["+vertexIJ(S.get(ss))[0]+","+vertexIJ(S.get(ss))[1]+","+"],");
                }
                System.out.println("]");
                System.out.print("Sink T: [");
                for(int ss=0;ss<T.size();ss++){
                    System.out.print("["+vertexIJ(T.get(ss))[0]+","+vertexIJ(T.get(ss))[1]+","+"],");
                }
                System.out.println("]");

                System.out.println("Image resolution: "+image.getWidth()+" x "+image.getHeight());

                long startT = System.currentTimeMillis();

                finalCuts = CountMinCuts.countCuts(G,S,T,0);

                System.out.println("Total Time taken: "+(System.currentTimeMillis()-startT)/1000);

                for(int xi=0;xi<image.getWidth();xi++){
                    for(int yi=0;yi<image.getHeight();yi++){
                        image.setRGB(xi,yi,imageOG.getRGB(xi,yi));
                    }
                }


                int[][] finalCut = finalCuts.get(0);
                for(int i=0;i<finalCut.length;i++){
                    int[] p1 = vertexIJ(finalCut[i][0]);
                    int[] p2 = vertexIJ(finalCut[i][1]);
                    image.setRGB(p1[0],p1[1],255);
                    image.setRGB(p2[0],p2[1],255);
                }
                if(finalCuts.size()>1) {
                    finalCut = finalCuts.get(1);
                    for (int i = 0; i < finalCut.length; i++) {
                        int[] p1 = vertexIJ(finalCut[i][0]);
                        int[] p2 = vertexIJ(finalCut[i][1]);
                        image.setRGB(p1[0], p1[1], 16776960);
                        image.setRGB(p2[0], p2[1], 16776960);
                    }
                }
                System.out.println("------------------x-------------------------x-----------------------------");

            }

            else if(e.getActionCommand().equalsIgnoreCase("Prob_Sample")){

                type = 1;

                menu2.setEnabled(true);
                menu.setEnabled(false);

                count++;
                scrib_s = false;
                scrib_t = false;

//                S.add(vertexNo(25,19));
//                T.add(vertexNo(17,19));

                System.out.print("Source S: [");
                for(int ss=0;ss<S.size();ss++){
                    System.out.print("["+vertexIJ(S.get(ss))[0]+","+vertexIJ(S.get(ss))[1]+","+"],");
                }
                System.out.println("]");
                System.out.print("Sink T: [");
                for(int ss=0;ss<T.size();ss++){
                    System.out.print("["+vertexIJ(T.get(ss))[0]+","+vertexIJ(T.get(ss))[1]+","+"],");
                }
                System.out.println("]");

                System.out.println("Image resolution: "+image.getWidth()+" x "+image.getHeight());

                String samplesString= JOptionPane.showInputDialog("No. of samples:");

                samples = Integer.parseInt(samplesString);

                long startT = System.currentTimeMillis();

                finalCuts = CountMinCuts.countCuts(G,S,T,samples);

                System.out.println("Total Time taken: "+(System.currentTimeMillis()-startT)/1000);



                for(int xi=0;xi<image.getWidth();xi++){
                    for(int yi=0;yi<image.getHeight();yi++){
                        image.setRGB(xi,yi,imageOG.getRGB(xi,yi));
                    }
                }

                pixCounts = new int[image.getWidth()][image.getHeight()];
                for(int i=0; i<finalCuts.size();i++){
                    int[][] finalCut = finalCuts.get(i);
                    for(int j=0;j<finalCut.length;j++) {
                        int[] p1 = vertexIJ(finalCut[j][0]);
                        int[] p2 = vertexIJ(finalCut[j][1]);
                        pixCounts[p1[0]][p1[1]]++;
                        pixCounts[p2[0]][p2[1]]++;
                    }
                }

                for(int xi=0;xi<image.getWidth();xi++){
                    for(int yi=0;yi<image.getHeight();yi++) {
                        if(pixCounts[xi][yi] != 0){
                            int value = (int)(((double)pixCounts[xi][yi]/samples)*256);
                            image.setRGB(xi,yi,16776960+(256-value));
                        }
                    }
                }


                System.out.println("------------------x-------------------------x-----------------------------");

            }
        }
    }


    void refreshPointer(MouseEvent e) {
        mouse_x = e.getX();
        mouse_y = e.getY();
        this.repaint();
    }

    public int[] RGB(int pixel) {
        int alpha = (pixel >> 24) & 0xff;
        int red = (pixel >> 16) & 0xff;
        int green = (pixel >> 8) & 0xff;
        int blue = (pixel) & 0xff;
        return new int[]{red,green,blue};
    }

//    public int RGBAValue(int red, int green, int blue, int alpha){
//        int r = red & 0xFF;
//        int g = green & 0xFF;
//        int b = blue & 0xFF;
//        int a = alpha & 0xFF;
//
//        return (r << 24) + (g << 16) + (b << 8) + (a);
//    }

    public int vertexNo(int i,int j) {
        return ((image.getHeight() * i) + j);
    }

    public int[] vertexIJ(int N) {
        return new int[]{N / image.getHeight(), N % image.getHeight()};
    }

    public double luminence(int[] rgb) {
//        return ((double)rgb[0]+rgb[1]+rgb[2])/3; // 1
//        return ((double)Math.max(Math.max(rgb[0],rgb[1]),rgb[2])+Math.min(Math.min(rgb[0],rgb[1]),rgb[2]))/2; // 2
//        return Math.max(Math.max(rgb[0],rgb[1]),rgb[2]); // 3
//        return Math.min(Math.min(rgb[0],rgb[1]),rgb[2]); // 4
//        int numShades= 100;
//        return ((int)(((double)rgb[0]+rgb[1]+rgb[2])/3 / (255 / ((double)numShades - 1))) + 0.5) * (255 / ((double)numShades - 1)); // 5
        return (0.2126 * rgb[0]) + (0.7152 * rgb[1]) + (0.0722 * rgb[2]); // 6
    }

    public double intensity(double diff) {
//        return Math.exp(-(diff))*Math.pow(10,22);
//        return Math.exp(-(Math.pow(diff,2)));

//        return (256-diff); // 3
//        return Math.pow((255-diff),8)+1; // 4
//        return 1/Math.pow(diff+1,2); // 5
//        return (Math.pow(((1 / (diff + 1)) * 1000), 4)* image.getWidth() * image.getHeight()); // 6
        return (Math.pow(((1 / (diff + 1)) * 1000), 8)* Math.pow(image.getWidth(),10) * Math.pow(image.getHeight(),10)); // 7
//        return Math.exp(-1*Math.sqrt(Math.pow(diff,3)));
//        return Math.exp(-1*diff);
//        return Math.exp(-(Math.pow(Math.sqrt(),2)/Math.pow(1,2)))
//        return Math.exp(-(Math.pow(Math.sqrt(diff),2)));
    }

    public int RGBInt(int avg){
        return (1<<24) | (avg<<16) | (avg<<8) | avg;
    }

    public ArrayList<Integer[]> findLine(int x0, int y0, int x1, int y1)
    {
        ArrayList<Integer[]> line = new ArrayList<>();

        int dx = Math.abs(x1 - x0);
        int dy = Math.abs(y1 - y0);

        int sx = x0 < x1 ? 1 : -1;
        int sy = y0 < y1 ? 1 : -1;

        int err = dx-dy;
        int e2;

        while (true)
        {
            line.add(new Integer[]{x0,y0});

            if (x0 == x1 && y0 == y1)
                break;

            e2 = 2 * err;
            if (e2 > -dy)
            {
                err = err - dy;
                x0 = x0 + sx;
            }

            else if (e2 < dx)
            {
                err = err + dx;
                y0 = y0 + sy;
            }
        }
        return line;
    }

    @Override
    public void mouseClicked(MouseEvent e) {

        if(count==1){

            if(x_prev == -1 && y_prev == -1){
                System.out.println("Initializing image");
                int maxDiff = 20;
                G = new adj_list(image.getHeight()*image.getWidth());
                int[][] grayScale = new int[image.getWidth()][image.getHeight()];

                for(int i=0; i<image.getWidth(); i++) {
                    for (int j = 0; j < image.getHeight(); j++) {
                        int nbrs[];
                        double weights[];
                        int rgb_ =image.getRGB(i,j);
                        int[] rgb__ = RGB(rgb_);
                        double lumz = luminence(rgb__);

                        grayScale[i][j] = RGBInt((int)lumz);

                        if((i==0 && j==0) || (i==image.getWidth()-1 && j==image.getHeight()-1) || (i==0 && j==image.getHeight()-1) || (i==image.getWidth()-1 && j==0)){
                            nbrs= new int[2];
                            weights = new double[2];
                        }
                        else if(i==0 || j==0 || i==image.getWidth()-1 || j==image.getHeight()-1){
                            nbrs = new int[3];
                            weights = new double[3];
                        }
                        else{
                            nbrs = new int[4];
                            weights = new double[4];
                        }
                        int[] pxl1 = RGB(image.getRGB(i, j));

                        int nbrI = 0;
                        if(j != 0){
                            int[] pxl2 = RGB(image.getRGB(i, j-1));
                            double diff = Math.abs(luminence(pxl1) - luminence(pxl2));
                            nbrs[nbrI] = vertexNo(i, j - 1);
                            if(diff>maxDiff){
                                diff=maxDiff;
                                weights[nbrI] = intensity(diff);
                            }
                            else{
                                weights[nbrI] = intensity(diff);
                            }
                            nbrI++;
                        }
                        if(i != 0){
                            int[] pxl2 = RGB(image.getRGB(i-1, j));
                            double diff = Math.abs(luminence(pxl1) - luminence(pxl2));
                            nbrs[nbrI] = vertexNo(i-1, j );
                            if(diff>maxDiff){
                                diff=maxDiff;
                                weights[nbrI] = intensity(diff);
                            }
                            else{
                                weights[nbrI] = intensity(diff);
                            }
                            nbrI++;
                        }

                        if(j != image.getHeight()-1){
                            int[] pxl2 = RGB(image.getRGB(i, j+1));
                            double diff = Math.abs(luminence(pxl1) - luminence(pxl2));
                            nbrs[nbrI] = vertexNo(i, j + 1);
                            if(diff>maxDiff){
                                diff=maxDiff;
                                weights[nbrI] = intensity(diff);
                            }
                            else{
                                weights[nbrI] = intensity(diff);
                            }
                            nbrI++;
                        }
                        if(i != image.getWidth()-1){
                            int[] pxl2 = RGB(image.getRGB(i+1, j));
                            double diff = Math.abs(luminence(pxl1) - luminence(pxl2));
                            nbrs[nbrI] = vertexNo(i+1, j);
                            if(diff>maxDiff){
                                diff=maxDiff;
                                weights[nbrI] = intensity(diff);
                            }
                            else{
                                weights[nbrI] = intensity(diff);
                            }
                            nbrI++;
                        }
                        G.insertVertex(vertexNo(i, j), nbrs, weights, false, new int[0][0]);
                    }
                }
                System.out.println("No. of vertices: "+G.num);


//                for(int i =0;i<grayScale.length;i++){
//                    for(int j=0;j<grayScale[0].length;j++){
//                        image.setRGB(i,j,grayScale[i][j]);
//                    }
//                }

                x_prev = e.getX();
                y_prev = e.getY();
                S.add(vertexNo(x_prev,y_prev));
                image.setRGB(x_prev,y_prev,65280);
                System.out.println(x+"--"+y);

            }
            else{
                x = e.getX();
                y = e.getY();
                System.out.println(x+"--"+y);
                ArrayList<Integer[]> line = findLine(x_prev,y_prev,x,y);
                int startI = -1;
                int sInd = -1;
                for(int i=0;i<S.size();i++) {
                    Integer xy = S.get(i);
                    sInd = -1;
                    int xy__2[] = vertexIJ(xy);
                    for(int j=0;j<line.size();j++){
                        Integer xy__1[] = line.get(j);
                        if(xy__1[0]==xy__2[0] && xy__1[1]==xy__2[1]){
                            sInd = j;
                        }
                    }
                    if(sInd>-1) {
                        startI = i;
                        break;
                    }
                }

                if(startI>-1){
                    S.subList(startI, S.size()).clear();
                }
                if(sInd==-1)
                    sInd=1;
                for(int i=sInd;i<line.size();i++){
                    Integer[] xy = line.get(i);
                    if(!S.contains(vertexNo(xy[0],xy[1]))){
                        S.add(vertexNo(xy[0],xy[1]));
                        image.setRGB(xy[0],xy[1],65280);
                    }
                    else{
                        System.out.println("---x----x----x----x----x---x---x---x----x---x----x----x---x--x----x---x----");
                    }
                }

                x_prev = x;
                y_prev = y;
            }

        }
        if(count==2){
            if(x_prev == -1 && y_prev == -1){
                x_prev = e.getX();
                y_prev = e.getY();
                T.add(vertexNo(x_prev,y_prev));
                image.setRGB(x_prev,y_prev,16711680);
                System.out.println(x_prev+"--"+y_prev);
            }
            else{

                x = e.getX();
                y = e.getY();
                System.out.println(x+"--"+y);
                ArrayList<Integer[]> line = findLine(x_prev,y_prev,x,y);
                int startI = -1;
                int sInd = -1;
                for(int i=0;i<T.size();i++) {
                    Integer xy = T.get(i);
                    sInd = -1;
                    int xy__2[] = vertexIJ(xy);
                    for(int j=0;j<line.size();j++){
                        Integer xy__1[] = line.get(j);
                        if(xy__1[0]==xy__2[0] && xy__1[1]==xy__2[1]){
                            sInd = j;
                        }
                    }
                    if(sInd>-1) {
                        startI = i;
                        break;
                    }
                }

                if(startI>-1){
                    T.subList(startI, T.size()).clear();
                }
                if(sInd==-1)
                    sInd=1;
                for(int i=sInd;i<line.size();i++){
                    Integer[] xy = line.get(i);
                    if(!T.contains(vertexNo(xy[0],xy[1]))){
                        T.add(vertexNo(xy[0],xy[1]));
                        image.setRGB(xy[0],xy[1],16711680);
                    }
                    else{
                        System.out.println("---x----x----x----x----x---x---x---x----x---x----x----x---x--x----x---x----");
                    }
                }
                x_prev = x;
                y_prev = y;
            }
        }
    }

    @Override
    public void mouseEntered(MouseEvent e) {
        // TODO Auto-generated method stub

    }

    @Override
    public void mouseExited(MouseEvent e) {
        // TODO Auto-generated method stub

    }

    @Override
    public void mousePressed(MouseEvent e) {
        // TODO Auto-generated method stub
    }

    @Override
    public void mouseReleased(MouseEvent e) {
        // TODO Auto-generated method stub
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        // TODO Auto-generated method stub
    }

    @Override
    public void mouseMoved(MouseEvent e) {
        refreshPointer(e);
    }


}

