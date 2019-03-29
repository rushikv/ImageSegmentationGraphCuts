/******************************************************************************
 *  Compilation:  javac CountMinCuts.java
 *  Execution:    java CountMinCuts
 *  Dependencies: IndexMinPQ.java
 *
 *  Count min cuts using Rachel Silva's algorithm
 *
 * Author: Rushik Vartak (rv9981@g.rit.edu)
 *
 ******************************************************************************/

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

class Node
{
    int id;
    double weight;
    Node ptr = null;
    Node left = null;
    Node right = null;
    Node leftFace = null;
    Node rightFace = null;
    ArrayList<Node> edges = new ArrayList<Node>();
    int edgeNo = 0;
    int[] edgeVertices;

    Node(int id_, double weight_){
        id = id_;
        weight = weight_;
    }

    String str()
    {
        if(ptr != null){
            return ""+id+" ("+weight+")->"+ptr.id;
        }
        else{
            return ""+id+" ("+weight+")->Null";
        }
    }
}


class adj_list
{

    ArrayList<ArrayList<Node>> list;
    double sum;
    int num;
    boolean visited[];
    int numEdges;
    int firstSV;
    int firstTV;


    adj_list(int num_){
        num = num_;
        list = new ArrayList<ArrayList<Node>>();
        for(int i=0;i<num;i++){
            list.add(new ArrayList<Node>());
        }
        sum = 0;
        visited= new boolean[num];
        numEdges = 0;
        firstSV = 0;
        firstTV = 0;
    }

    void insertVertex(int V, int neighbors[], double weights[], boolean dual, int edgeVs[][])
    {
        int count = 0;
        visited[V] = true;

        for(int i=0; i< neighbors.length; i++){
            int nbr = neighbors[i];
            double weight = weights[i];
            Node ownNode = null;
            Node nbrNode = null;
            if(!visited[nbr]){
                numEdges++;
                nbrNode = new Node(nbr,weight);
                ownNode = new Node(V, weight);
                nbrNode.ptr = ownNode;
                ownNode.ptr = nbrNode;
                nbrNode.edgeNo = count;
                list.get(V).add(nbrNode);
                sum=sum+weight;
            }
            else{
                if(dual == true){
                    ownNode = getUniqueNodeFromList(V,edgeVs[i],list.get(nbr));
                }
                else{
                    ownNode = getNodeFromList(V,list.get(nbr));
                }
                if(ownNode==null) {
                    System.out.println(V+"--("+edgeVs[i][0]+","+edgeVs[i][1]+")---"+nbr);
                }
                nbrNode = ownNode.ptr;
                nbrNode.edgeNo = count;
                list.get(V).add(nbrNode);
            }

            if(dual == true){
                nbrNode.edgeVertices = edgeVs[i];
                int revNbrs[] = {edgeVs[i][1],edgeVs[i][0]};
                ownNode.edgeVertices = revNbrs;
                count+=1;
            }
        }
    }

    Node getUniqueNodeFromList(int V, int[] edgeVs, ArrayList<Node> list){
        for(int i=0; i < list.size(); i++ ){
            Node curr = list.get(i);
            if(curr.id == V){
                if((curr.edgeVertices[0] == edgeVs[0] && curr.edgeVertices[1] == edgeVs[1]) || (curr.edgeVertices[0] == edgeVs[1] && curr.edgeVertices[1] == edgeVs[0])) {
                    return curr;
                }
            }
        }
        return null;
    }

    Node getNodeFromList(int V, ArrayList<Node> list){
        for(int i=0; i < list.size(); i++ ){
            Node curr = list.get(i);
            if(curr.id == V){
                return curr;
            }
        }
        return null;
    }

    int getIndexFromList(int V, ArrayList<Node> list){
        for(int i=0; i < list.size(); i++ ){
            if(list.get(i).id == V){
                return i;
            }
        }
        return -1;
    }

    void copy (adj_list G_copy){
        for(int key=0; key<num; key++){
            int[] nbr_ids = new int[list.get(key).size()];
            double[] nbr_weights = new double[list.get(key).size()];

            for(int i=0; i<list.get(key).size();i++){
//                G_copy.list.get(key).add(list.get(key).get(i));
                nbr_ids[i] = list.get(key).get(i).id;
                nbr_weights[i] = list.get(key).get(i).weight;
            }
            G_copy.insertVertex(key, nbr_ids, nbr_weights, false, new int[0][0]);
        }
//        G_copy.sum = sum;
//        G_copy.numEdges = numEdges;
    }

    ArrayList<Integer> createFace(ArrayList<Integer> SList, int vertex){
        ArrayList<Integer> newVertices = new ArrayList<Integer>();
        ArrayList<Node> newVs = new ArrayList<Node>();
        int ctr = 0;
        int S_next, S_prev, i_next, i_prev;

        if(SList.size()==1){
            int S = SList.get(0);
            for (int i = 0; i < list.get(S).size(); i++) {
                Node nbr = list.get(S).get(i);
                if(!SList.contains(nbr.id)) {
                    Node ownNode = nbr.ptr;
                    newVs.add(ownNode);
                    newVertices.add(num + ctr);
                    ctr += 1;
                }
            }
        }
        else{
            int S = SList.get(0);

            S_next = SList.get(1);
            i_next = getIndexFromList(S_next, list.get(S));

            for (int i = i_next; i < list.get(S).size(); i++) {
                Node nbr = list.get(S).get(i);
                if(!SList.contains(nbr.id)) {
                    Node ownNode = nbr.ptr;
                    newVs.add(ownNode);
                    newVertices.add(num + ctr);
                    ctr += 1;
                }
            }
            for (int i = 0; i < i_next; i++) {
                Node nbr = list.get(S).get(i);
                if(!SList.contains(nbr.id)) {
                    Node ownNode = nbr.ptr;
                    newVs.add(ownNode);
                    newVertices.add(num + ctr);
                    ctr += 1;
                }
            }

            for(int z=1; z<SList.size()-1; z++) {
                S = SList.get(z);
                S_prev = SList.get(z-1);
                S_next = SList.get(z+1);
                i_prev = getIndexFromList(S_prev, list.get(S));
                i_next = getIndexFromList(S_next, list.get(S));

                if(i_prev<i_next){
                    for (int i = i_prev; i < i_next; i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                }
                else{
                    for (int i = i_prev; i < list.get(S).size(); i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                    for (int i = 0; i < i_next; i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                }
            }

            S = SList.get(SList.size()-1);
            S_prev = SList.get(SList.size()-2);
            i_prev = getIndexFromList(S_prev, list.get(S));
            for (int i = i_prev; i < list.get(S).size(); i++) {
                Node nbr = list.get(S).get(i);
                if(!SList.contains(nbr.id)) {
                    Node ownNode = nbr.ptr;
                    newVs.add(ownNode);
                    newVertices.add(num + ctr);
                    ctr += 1;
                }
            }
            for (int i = 0; i < i_prev; i++) {
                Node nbr = list.get(S).get(i);
                if(!SList.contains(nbr.id)) {
                    Node ownNode = nbr.ptr;
                    newVs.add(ownNode);
                    newVertices.add(num + ctr);
                    ctr += 1;
                }
            }

            for(int z=SList.size()-2; z>0; z--) {
                S = SList.get(z);
                S_prev = SList.get(z-1);
                S_next = SList.get(z+1);
                i_prev = getIndexFromList(S_prev, list.get(S));
                i_next = getIndexFromList(S_next, list.get(S));

                if(i_next<i_prev){
                    for (int i = i_next; i < i_prev; i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                }
                else{
                    for (int i = i_next; i < list.get(S).size(); i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                    for (int i = 0; i < i_prev; i++) {
                        Node nbr = list.get(S).get(i);
                        if(!SList.contains(nbr.id)) {
                            Node ownNode = nbr.ptr;
                            newVs.add(ownNode);
                            newVertices.add(num + ctr);
                            ctr += 1;
                        }
                    }
                }
            }
        }

        if(vertex == 0) {
            firstSV = num;
        }
        else if(vertex == 1) {
            firstTV = num;
        }

        for(int i=0; i<newVs.size(); i ++){
            list.add(new ArrayList<Node>());
        }

        ctr = 1;
        Node firstV = newVs.get(0);
        firstV.id = newVertices.get(0);
        list.get(firstV.id).add(firstV.ptr);
        Node nextNode = new Node((num + (ctr%newVs.size())), sum);
        Node ownNode = new Node(firstV.id, sum);
        nextNode.ptr = ownNode;
        ownNode.ptr = nextNode;
        list.get(firstV.id).add(nextNode);
        Node previous = ownNode;
        numEdges+=1;

        for(int i = 1; i<newVs.size();i++){
            ctr+=1;
            Node curr=newVs.get(i);
            curr.id = newVertices.get(i);
            list.get(curr.id).add(curr.ptr);
            nextNode = new Node(num + (ctr%newVs.size()), sum);
            ownNode = new Node(newVs.get(i).id, sum);
            nextNode.ptr = ownNode;
            ownNode.ptr = nextNode;
            list.get(curr.id).add(nextNode);
            list.get(curr.id).add(previous);
            previous = ownNode;
            numEdges += 1;
        }

        list.get(firstV.id).add(previous);
        for(int z=0; z<SList.size(); z++) {
            int S = SList.get(z);
            list.set(S, new ArrayList<Node>());
        }
        num+=newVs.size();

        return newVertices;
    }


    void linkNodes(){
        for(int currVKey=0; currVKey<num; currVKey++){
            ArrayList<Node> Es = list.get(currVKey);
            for(int i=0; i<Es.size(); i++){
                if(i==0){
                    Es.get(i).left = Es.get(Es.size()-1);
                }
                else{
                    Es.get(i).left = Es.get(i-1);
                }
                Es.get(i).right = Es.get((i+1)%Es.size());
                Es.get(i).edgeNo = i;
            }
        }
    }

    void dual(adj_list Dual){
        int count = 0;
        ArrayList<Node> faceList = new ArrayList<Node>();

        for(int currVKey = 0; currVKey<num; currVKey++){
            ArrayList<Node> Es = list.get(currVKey);
            for (Node currEdge:Es) {
                Node currNode = currEdge;
                boolean broke = false;
                ArrayList<Node> faceEdges = new ArrayList<Node>();
                faceEdges.add(currNode);
                int s = 0;
                int t = 0;

                if(currNode.id >= firstSV && currNode.id < firstTV){
                    s++;
                }
                else if (currNode.id >= firstTV){
                    t++;
                }

                while (currNode.id != currVKey){
                    if(currNode.leftFace == null){
                        currNode = currNode.ptr.right;
                        faceEdges.add(currNode);
                        if(currNode.id >= firstSV && currNode.id < firstTV) {
                            s++;
                        }
                        else if(currNode.id >= firstTV) {
                            t++;
                        }
                    }
                    else{
                        broke = true;
                        break;
                    }
                }

                if(!broke){
                    if(s==faceEdges.size()){
                        Dual.firstSV = count;
                    }
                    else if(t==faceEdges.size()){
                        Dual.firstTV = count;
                    }
                    Node face = new Node(count, 0);
                    count++;
                    for(int i=faceEdges.size()-1;i>=0;i--){
                        Node node = faceEdges.get(i);
                        face.edges.add(node);
                        node.leftFace = face;
                        node.ptr.rightFace = face;
                    }
                    faceList.add(face);
                }

            }
        }

        for (Node face: faceList) {
            int[] nbrs = new int[face.edges.size()];
            double[] weights = new double[face.edges.size()];
            int[][] edgeVs = new int[face.edges.size()][2];

            for(int i=0; i < face.edges.size(); i++){
                nbrs[i] = face.edges.get(i).rightFace.id;
                weights[i] = face.edges.get(i).weight;
                edgeVs[i][0] = face.edges.get(i).id;
                edgeVs[i][1] = face.edges.get((i+1)%face.edges.size()).id;
            }
            Dual.insertVertex(face.id, nbrs, weights, true, edgeVs);
        }
    }
}


class DijRetObj{
    double[] visited;
    double[] count;
    ArrayList<ArrayList<Integer[]>> path;
    int start;
    int end;
    DijRetObj(double[] visited_, double[] count_, ArrayList<ArrayList<Integer[]>> path_, int start_, int end_){
        visited = visited_;
        count = count_;
        path = path_;
        start = start_;
        end = end_;
    }
}

public class CountMinCuts {

    private static DijRetObj DijkstraCount(adj_list graph, int initial, int dest){
        double[] cost = new double[graph.num];
        double[] count = new double[graph.num];
        IndexMinPQ<Double> queue = new IndexMinPQ<>(graph.num);
        queue.insert(initial,0.0);
        ArrayList<ArrayList<Integer[]>> path = new ArrayList<ArrayList<Integer[]>>();
        for(int i=0;i<graph.num;i++){
            path.add(new ArrayList<Integer[]>());
        }
        int nodesCount = 0;
        count[initial] = 1;

        while(nodesCount!=graph.num && !queue.isEmpty()){
            double current_weight = queue.minKey();
            int min_node = queue.delMin();

            nodesCount++;

            for (Node v: graph.list.get(min_node)) {
                double weight = current_weight + v.weight;
                if(weight==cost[v.id]){
                    Integer[] arr = {min_node, v.edgeNo};
                    path.get(v.id).add(arr);
                    count[v.id]+=count[min_node];
                }
                else if(cost[v.id] == 0 || weight < cost[v.id]) {
                    cost[v.id] = weight;

                    if(!queue.contains(v.id)){
                        queue.insert(v.id, weight);
                    }
                    else{
                        queue.decreaseKey(v.id, weight);
                    }
                    Integer[] arr = {min_node, v.edgeNo};
                    path.set(v.id,new ArrayList<Integer[]>());
                    path.get(v.id).add(arr);
                    count[v.id] = count[min_node];
                }

            }
        }
        return new DijRetObj(cost, count, path, initial, dest);
    }

    public static ArrayList<Integer[]> RMSP(adj_list Dual, ArrayList<ArrayList<Integer[]>> DijTree, int source, int dest){
        Integer[] IntArr_2;
        ArrayList<Integer[]> RMPath = new ArrayList<Integer[]>();
        Integer[] IntArr = {source,Dual.list.get(DijTree.get(source).get(0)[0]).get(DijTree.get(source).get(0)[1]).ptr.edgeNo};
        RMPath.add(IntArr);

        int curr = DijTree.get(source).get(0)[0];
        int prevEdge = DijTree.get(source).get(0)[1];
        Integer[] currEdge;
        while(curr!=dest){
            currEdge = DijTree.get(curr).get(0);
            for(int i=1;i<DijTree.get(curr).size();i++){
                int currEdgeNo = Dual.list.get(currEdge[0]).get(currEdge[1]).ptr.edgeNo;
                    if (Dual.list.get(DijTree.get(curr).get(i)[0]).get(DijTree.get(curr).get(i)[1]).ptr.edgeNo > currEdgeNo && Dual.list.get(DijTree.get(curr).get(i)[0]).get(DijTree.get(curr).get(i)[1]).ptr.edgeNo < prevEdge) {
                        currEdge = DijTree.get(curr).get(i);
                    }

            }
            IntArr_2 = new Integer[]{curr,Dual.list.get(currEdge[0]).get(currEdge[1]).ptr.edgeNo};
            RMPath.add(IntArr_2);
            curr = currEdge[0];
            prevEdge = currEdge[1];
        }
        return RMPath;
    }

    public static ArrayList<Integer[]> LMSP(adj_list Dual, ArrayList<ArrayList<Integer[]>> DijTree, int source, int dest){
        Integer[] IntArr_2;
        ArrayList<Integer[]> LMPath = new ArrayList<Integer[]>();
        Integer[] IntArr = {source,Dual.list.get(DijTree.get(source).get(0)[0]).get(DijTree.get(source).get(0)[1]).ptr.edgeNo};
        LMPath.add(IntArr);

        int curr = DijTree.get(source).get(0)[0];
        int prevEdge = DijTree.get(source).get(0)[1];
        Integer[] currEdge;
        while(curr!=dest){
            currEdge = DijTree.get(curr).get(0);
            for(int i=1;i<DijTree.get(curr).size();i++){
                int currEdgeNo = Dual.list.get(currEdge[0]).get(currEdge[1]).ptr.edgeNo;
                    if (Dual.list.get(DijTree.get(curr).get(i)[0]).get(DijTree.get(curr).get(i)[1]).ptr.edgeNo < currEdgeNo && Dual.list.get(DijTree.get(curr).get(i)[0]).get(DijTree.get(curr).get(i)[1]).ptr.edgeNo > prevEdge) {
                        currEdge = DijTree.get(curr).get(i);
                    }
            }
            IntArr_2 = new Integer[]{curr,Dual.list.get(currEdge[0]).get(currEdge[1]).ptr.edgeNo};
            LMPath.add(IntArr_2);
            curr = currEdge[0];
            prevEdge = currEdge[1];
        }

        return LMPath;
    }

    public static int[] transformRMSP(adj_list Dual, ArrayList<Integer[]> RMPath){
        for(int i=0;i<RMPath.size()-1;i++){
            Dual.list.add(new ArrayList<Node>());
        }

        int[] newNodes = new int[RMPath.size()-1];
        int count = 0;

        for(int i=1; i<RMPath.size(); i++){
            Integer[] prev = RMPath.get(i-1);
            Integer[] curr = RMPath.get(i);
            Node V = Dual.list.get(prev[0]).get(prev[1]).ptr;
            Node nextV = Dual.list.get(curr[0]).get(curr[1]);
            Node currV = V.left;
            newNodes[i-1]=Dual.num+count;

            while(currV.id != nextV.id && currV.edgeNo != nextV.edgeNo){
                currV.ptr.id = Dual.num + count;
                Dual.list.get(currV.ptr.id).add(0,currV);
                currV.left.right = currV.right;
                currV.right.left = currV.left;
                if(currV.edgeNo < curr[1]){
                    RMPath.get(i)[1]--;
                    curr[1]=RMPath.get(i)[1];
                }
                for(int ind=currV.edgeNo+1; ind<Dual.list.get(curr[0]).size(); ind++) {
                    Dual.list.get(curr[0]).get(ind).edgeNo = ind - 1;
                }

                Dual.list.get(curr[0]).remove(currV.edgeNo);
                currV = currV.left;
            }
            count++;
        }

        Dual.num += newNodes.length;

        return newNodes;
    }

    public static Integer[] probRandomList(ArrayList<Integer[]> list, double[] count, double sum){

        double rand = ThreadLocalRandom.current().nextDouble(1, sum);
        double tempSum = 0;
        for(int i=0;i<list.size();i++){
            tempSum+=count[i];
            if(tempSum>=rand){
                return list.get(i);
            }
        }
        return null;
    }

    public static DijRetObj probRandomList(ArrayList<DijRetObj> list, ArrayList<Double> count, double sum){

        double rand = ThreadLocalRandom.current().nextDouble(1, sum);
        double tempSum = 0;
        for(int i=0;i<list.size();i++){
            tempSum+=count.get(i);
            if(tempSum>=rand){
                return list.get(i);
            }
        }
        return null;
    }

    public static ArrayList<Integer[]> probPath(ArrayList<ArrayList<Integer[]>> path, double[] count, int initial, int dest){
        int curr = initial;
        ArrayList<Integer[]> chosenPath = new ArrayList<Integer[]>();

        while(curr!=dest){
            if(path.get(curr).size() == 1){
                chosenPath.add(path.get(curr).get(0));
                curr = path.get(curr).get(0)[0];
            }
            else{
                double countList[] = new double[path.get(curr).size()];
                for(int i=0; i<path.get(curr).size(); i++){
                    countList[i] = count[path.get(curr).get(i)[0]];
                }
                Integer[] probChosen = probRandomList(path.get(curr),countList,count[curr]);
                chosenPath.add(probChosen);
                curr = probChosen[0];
            }
        }
        return chosenPath;
    }


    public static ArrayList<int[][]> countCuts(adj_list G, ArrayList<Integer> S, ArrayList<Integer> T, int samples){
        adj_list G_1 = new adj_list(G.num);

        G.copy(G_1);
        System.out.println("------------------G_1 created-------------------");
        ArrayList<Integer> sourceVs = G_1.createFace(S,0);
        System.out.println("------------------S face made-------------------");
        ArrayList<Integer> sinkVs = G_1.createFace(T,1);
        System.out.println("------------------T face made-------------------");
        G_1.linkNodes();

        adj_list Dual = new adj_list(G_1.numEdges - G_1.num + 4);
        G_1.dual(Dual);
        Dual.linkNodes();
        System.out.println("------------------Dual created-------------------");

        DijRetObj dijResult = DijkstraCount(Dual, Dual.firstTV, Dual.firstSV);
        ArrayList<Integer[]> rightMostPath = RMSP(Dual, dijResult.path, Dual.firstSV, Dual.firstTV);
        System.out.println("Right Most shortest path:");
        for(int i=0;i<rightMostPath.size();i++){
            System.out.print("["+rightMostPath.get(i)[0]+","+rightMostPath.get(i)[1]+"], ");
        }
        System.out.println();
        int[] newNodes = transformRMSP(Dual,rightMostPath);
        Dual.linkNodes();



        double minDist = Double.MAX_VALUE;
        double minCount = 0;

        DijRetObj firstObj = null;
        DijRetObj lastObj = null;
        int firstStart =-1;
        int firstEnd = -1;
        int lastStart =-1;
        int lastEnd = -1;


        ArrayList<DijRetObj> shortestDijObjs = new ArrayList<DijRetObj>();
        ArrayList<Double> chosenCounts = new ArrayList<Double>();

        int perc[] = new int[]{0, 25,50,75,100};
        int percI = 0;


        System.out.println("Size of RMSP (d):"+ rightMostPath.size());

        for(int i=1; i<rightMostPath.size(); i++){
            DijRetObj currDijResult = DijkstraCount(Dual, newNodes[i-1], rightMostPath.get(i)[0]);
            double[] count = currDijResult.count;
            double[] dist = currDijResult.visited;
            if((i*100.0)/(double)rightMostPath.size()>=perc[percI]){
                System.out.println(perc[percI]+"% completed!");
                percI++;
            }

            if(count[rightMostPath.get(i)[0]]==0){
                continue;
            }

            if(dist[rightMostPath.get(i)[0]]==minDist){
                minCount += count[rightMostPath.get(i)[0]];

                if(samples==0) {
                    lastObj = currDijResult;
                    lastStart = i;
                    lastEnd = i - 1;
                }
                else{
                    shortestDijObjs.add(currDijResult);
                    chosenCounts.add(count[rightMostPath.get(i)[0]]);
                }
            }
            else if(dist[rightMostPath.get(i)[0]]<minDist){
                minCount = count[rightMostPath.get(i)[0]];
                minDist = dist[rightMostPath.get(i)[0]];
                if(samples==0) {
                    firstObj = currDijResult;
                    lastObj = currDijResult;
                    firstStart = i;
                    firstEnd = i - 1;
                    lastStart = i;
                    lastEnd = i - 1;
                }
                else{
                    shortestDijObjs = new ArrayList<DijRetObj>();
                    shortestDijObjs.add(currDijResult);
                    chosenCounts = new ArrayList<Double>();
                    chosenCounts.add(count[rightMostPath.get(i)[0]]);
                }
            }
        }
        System.out.println("100% completed!");

        System.out.println("Count: "+minCount);

        ArrayList<int[][]> returnResult = new ArrayList<int[][]>();
        ArrayList<Integer[]> finalPath;
        ArrayList<ArrayList<Integer[]>> finalPath_List;

        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++===");
        if(minCount>1) {
            if(samples==0) {
                finalPath_List = new ArrayList<ArrayList<Integer[]>>();
                finalPath_List.add(LMSP(Dual, firstObj.path, rightMostPath.get(firstStart)[0], newNodes[firstEnd]));
                finalPath_List.add(RMSP(Dual, lastObj.path, rightMostPath.get(lastStart)[0], newNodes[lastEnd]));
                int[][] planarPath_1 = new int[finalPath_List.get(0).size()][2];
                int[][] planarPath_2 = new int[finalPath_List.get(1).size()][2];


                for (int i = finalPath_List.get(0).size() - 1; i >= 0; i--) {
                    int[] E = Dual.list.get(finalPath_List.get(0).get(i)[0]).get(finalPath_List.get(0).get(i)[1]).edgeVertices;

                    if (E[0] >= sourceVs.get(0) && E[0] <= sourceVs.get(sourceVs.size() - 1)) {
                        E[0] = S.get(0);
                    } else if (E[0] >= sinkVs.get(0) && E[0] <= sinkVs.get(sinkVs.size() - 1)) {
                        E[0] = T.get(0);
                    }

                    if (E[1] >= sourceVs.get(0) && E[1] <= sourceVs.get(sourceVs.size() - 1)) {
                        E[1] = S.get(0);
                    } else if (E[1] >= sinkVs.get(0) && E[1] <= sinkVs.get(sinkVs.size() - 1)) {
                        E[1] = T.get(0);
                    }

                    planarPath_1[finalPath_List.get(0).size() - i - 1][0] = E[0];
                    planarPath_1[finalPath_List.get(0).size() - i - 1][1] = E[1];
                }

                for (int i = finalPath_List.get(1).size() - 1; i >= 0; i--) {
                    int[] E = Dual.list.get(finalPath_List.get(1).get(i)[0]).get(finalPath_List.get(1).get(i)[1]).edgeVertices;

                    if (E[0] >= sourceVs.get(0) && E[0] <= sourceVs.get(sourceVs.size() - 1)) {
                        E[0] = S.get(0);
                    } else if (E[0] >= sinkVs.get(0) && E[0] <= sinkVs.get(sinkVs.size() - 1)) {
                        E[0] = T.get(0);
                    }

                    if (E[1] >= sourceVs.get(0) && E[1] <= sourceVs.get(sourceVs.size() - 1)) {
                        E[1] = S.get(0);
                    } else if (E[1] >= sinkVs.get(0) && E[1] <= sinkVs.get(sinkVs.size() - 1)) {
                        E[1] = T.get(0);
                    }

                    planarPath_2[finalPath_List.get(1).size() - i - 1][0] = E[0];
                    planarPath_2[finalPath_List.get(1).size() - i - 1][1] = E[1];
                }

                System.out.print("Leftmost Final path: ");

                for (int i = 0; i < planarPath_1.length; i++) {
                    System.out.print("[" + planarPath_1[i][0] + "," + planarPath_1[i][1] + "],");
                }
                System.out.println();

                System.out.print("Rightmost Final path: ");

                for (int i = 0; i < planarPath_2.length; i++) {
                    System.out.print("[" + planarPath_2[i][0] + "," + planarPath_2[i][1] + "],");
                }
                System.out.println();

                returnResult.add(planarPath_1);
                returnResult.add(planarPath_2);

                return returnResult;
            }
            else{
                System.out.println("No. of Dij Obs stored: "+shortestDijObjs.size());
                for(int Zi=0;Zi<samples;Zi++){
                    DijRetObj probObj = probRandomList(shortestDijObjs,chosenCounts,minCount);
                    finalPath = probPath(probObj.path, probObj.count, probObj.end, probObj.start);
//                    System.out.print("Chosen path: ");
//
//                    for(int i=0; i <finalPath.size(); i ++){
//                        System.out.print("["+finalPath.get(i)[0]+","+finalPath.get(i)[1]+"],");
//                    }
//                    System.out.println();

                    int[][] planarPath = new int[finalPath.size()][2];

                    for(int i=finalPath.size()-1; i>=0; i--){
                        int[] E = Dual.list.get(finalPath.get(i)[0]).get(finalPath.get(i)[1]).edgeVertices;

                        if(E[0] >= sourceVs.get(0) && E[0] <= sourceVs.get(sourceVs.size()-1)){
                            E[0] = S.get(0);
                        }
                        else if(E[0] >= sinkVs.get(0) && E[0] <= sinkVs.get(sinkVs.size()-1)){
                            E[0] = T.get(0);
                        }

                        if(E[1] >= sourceVs.get(0) && E[1] <= sourceVs.get(sourceVs.size()-1)){
                            E[1] = S.get(0);
                        }
                        else if(E[1] >= sinkVs.get(0) && E[1] <= sinkVs.get(sinkVs.size()-1)){
                            E[1] = T.get(0);
                        }
                        planarPath[finalPath.size()-i-1][0]=E[0];
                        planarPath[finalPath.size()-i-1][1]=E[1];
                    }

//                    System.out.print("Final path: ");
//
//                    for(int i=0; i <planarPath.length; i ++){
//                        System.out.print("["+planarPath[i][0]+","+planarPath[i][1]+"],");
//                    }
//                    System.out.println();

                    returnResult.add(planarPath);

                }
                return returnResult;
            }

        }
        else{
            finalPath = probPath(firstObj.path, firstObj.count, rightMostPath.get(firstStart)[0], newNodes[firstEnd]);
            System.out.print("Chosen path: ");

            for(int i=0; i <finalPath.size(); i ++){
                System.out.print("["+finalPath.get(i)[0]+","+finalPath.get(i)[1]+"],");
            }
            System.out.println();

            int[][] planarPath = new int[finalPath.size()][2];

            for(int i=finalPath.size()-1; i>=0; i--){
                int[] E = Dual.list.get(finalPath.get(i)[0]).get(finalPath.get(i)[1]).edgeVertices;

                if(E[0] >= sourceVs.get(0) && E[0] <= sourceVs.get(sourceVs.size()-1)){
                    E[0] = S.get(0);
                }
                else if(E[0] >= sinkVs.get(0) && E[0] <= sinkVs.get(sinkVs.size()-1)){
                    E[0] = T.get(0);
                }

                if(E[1] >= sourceVs.get(0) && E[1] <= sourceVs.get(sourceVs.size()-1)){
                    E[1] = S.get(0);
                }
                else if(E[1] >= sinkVs.get(0) && E[1] <= sinkVs.get(sinkVs.size()-1)){
                    E[1] = T.get(0);
                }
                planarPath[finalPath.size()-i-1][0]=E[0];
                planarPath[finalPath.size()-i-1][1]=E[1];
            }

            System.out.print("Final path: ");

            for(int i=0; i <planarPath.length; i ++){
                System.out.print("["+planarPath[i][0]+","+planarPath[i][1]+"],");
            }
            System.out.println();

            returnResult.add(planarPath);
            return returnResult;
        }
    }

    public static void main(String args[]){
        adj_list G = new adj_list(10);
        G.insertVertex(0, new int[]{1,8,5}, new double[]{1.0, 2.0, 1.0}, false, new int[0][0]);
        G.insertVertex(1, new int[]{2,6,8,0}, new double[]{1.0, 1.0, 1.0, 1.0}, false, new int[0][0]);
        G.insertVertex(2, new int[]{3,9,7,1}, new double[]{2.0, 2.0, 1.0, 1.0}, false, new int[0][0]);
        G.insertVertex(3, new int[]{4,9,2}, new double[]{1.0, 3.0, 2.0}, false, new int[0][0]);
        G.insertVertex(4, new int[]{5,7,9,3}, new double[]{2.0, 3.0, 1.0, 1.0}, false, new int[0][0]);
        G.insertVertex(5, new int[]{0,8,6,4}, new double[]{1.0, 1.0, 3.0, 2.0}, false, new int[0][0]);
        G.insertVertex(6, new int[]{8,1,7,5}, new double[]{1.0, 1.0, 2.0, 3.0}, false, new int[0][0]);
        G.insertVertex(7, new int[]{6,2,9,4}, new double[]{2.0, 1.0, 1.0, 3.0}, false, new int[0][0]);
        G.insertVertex(8, new int[]{0,1,6,5}, new double[]{2.0, 1.0, 1.0, 1.0}, false, new int[0][0]);
        G.insertVertex(9, new int[]{7,2,3,4}, new double[]{1.0, 2.0, 3.0, 1.0}, false, new int[0][0]);

        ArrayList<Integer> S = new ArrayList<Integer>(){};
        S.add(8);
        ArrayList<Integer> T = new ArrayList<Integer>(){};
        T.add(9);
        countCuts(G,S,T,0);
    }
}
