Array.prototype.intersection = function() {
    var i, l = this.length, k = arguments, o = [];
    for(i=0; i<l;i+=1){
      if(k[0].includes(this[i]))
        o.push(this[i])
    }
    return o;
};

Array.prototype.getGeneString = function(){
  var str="",l = this.length, k = arguments;
  for(var i=0;i<l;i++){
    str += this[i]+k[0]
  }

  return str.slice(0,str.length-1);
}


function createGeneInfo(){
  gene_info = {};
  go_inf.forEach(function(d,i){
  var index = i;
  d["genes"].forEach(function(d,i){
    if(gene_info[d]!=undefined){
        var str =   ";"  +  index.toString();
        gene_info[d] = gene_info[d] + str;
    }
    else
       gene_info[d] = index.toString();
  });
});
}

function createInteractionGenes(go_inf){
   dic = {};
   var size = go_inf.length;
   for(var i = 0;i < size;i++){
      for(var j = 0;j < size;j++){
        if(i == j){
          dic[i+"-"+j] = ""
        }else{
          var genes_i = go_inf[i]["genes"];
          var genes_j = go_inf[j]["genes"];
          var intersectionGenes = genes_i.intersection(genes_j);
          

          dic[i+"-"+j] = intersectionGenes.getGeneString("|");
        }
      }
   }
}

function GOFilter(){
  var go_filter = [];
  go_inf.map(function(d,i){
      if(Number(d.pVal) < 0.05){
        go_filter.push(d);//remove
      }
  });
  go_inf = go_filter;
  nodesSize = go_inf.length;
};

function getMinPval(){
  var min = Number(go_inf[0].pVal);
  for(i=0;i<nodesSize;i++){
    if(Number(go_inf[i].pVal) < min){
      min = Number(go_inf[i].pVal);

    }
  }
  return min;
}

function parentGO(goid) {
  goids=[];
  var queue = new Queue();
  var go_level = [goid,1];
  goids.push(go_level);
  queue.enqueue(go_level);

  while(!queue.isEmpty()){
    go_level = queue.dequeue();

    go = go_level[0];
    level = go_level[1];

    current = GO(go);
    if(current!=undefined){
      var parents = current['p'];//there is big difference without var
      for (var i = 0, count = parents.length ; i < count; i ++ ) {
        if(parents[i]!=undefined)
        {
          goid = parents[i];
          go_level = [goid,level+1];

          goids.push(go_level);
          queue.enqueue(go_level);
        }
      }
    }
  }

  return goids;
}



//creates a D3 structure that can be inputted into the force-directed graph layout
function createD3Structure(goids, target) {

  goids = getUniqueGoid(goids);

  nodes = [];
  links = [];
  
  nodeIndex = {};
  levels = goids[goids.length-1][1]; //total level
  level_width ={};

  for (var i = 0, count = goids.length ; i < count ; i ++ ) {
    go_level = goids[i];
    goid = go_level[0];
    level = go_level[1];
    node = {};
    node['name'] = GO(goid)['n'];
    node['id'] = goid;
    node['r'] = 10;
    node['is'] = 'node';
    nodeIndex[goid] = i;
    

    //get the number of GO on each level
    if(level_width[level]==undefined){
      level_width[level]=1;
      node['position'] = 1;
    }
    else{
      level_width[level]++;
      node['position'] = level_width[level];
    }

    nodes.push(node);
  }

  for (var i = 0, count = goids.length ; i < count ; i ++ ) {
    go_level = goids[i];
    goid = go_level[0];
    level = go_level[1];

    parents = GO(goid)['p'];

    if (parents.length > 0 ) {
      for (var a = 0, countA = parents.length ; a < countA ; a ++ ) {
        parent = parents[a];
        pgoid = parent;
        pIndex = nodeIndex[pgoid];

     /*   position = nodes[pIndex]['position']; 
        nodes[pIndex].fixed = true;
        nodes[pIndex].x = width*(position/level_width[level]*0.5);
        nodes[pIndex].y = height*(1-(level)/levels);*/

        if (pgoid == target) {
          nodes[i]['is'] = 'child';
          nodes[i]['r'] = 30;
        }

        if (pIndex != undefined) {
          link = {};
          link['source'] = pIndex;
          link['target'] = i;
          link['value'] = 5;
          //link['type'] = prel;
          links.push(link)
        }
      }
    } else {        
      nodes[i].fixed = true; //the top node
      nodes[i].x = width / 2;
      nodes[i].y = 50;
      nodes[i]['is'] = 'root';
    }
  }


  nodes[nodeIndex[target]].fixed = true;
  nodes[nodeIndex[target]].x = width / 2;
  nodes[nodeIndex[target]].y = height - 20;
  nodes[nodeIndex[target]]['is'] = 'target';

  return {'nodes' : nodes, 'links' : links};
}      

function GO(goid) {
  return goNodes[goid];
}

function recursiveAncestor(goid) {
  goids = [];
  goids.push(goid);

  var data = GO(goid);

  if (data['p'].length > 0) {
    for (var i = 0, count = data['p'].length ; i < count ; i ++ ) {
      pgoid = data['p'][i];
      goids = goids.concat(recursiveAncestor(pgoid));
    }
  }

  return goids
}

/*    function immediateChildren(goid) {
  current = GO(goid);
  childrenIndeces = current['c'];

  goids = [];
  for (var i = 0, count = childrenIndeces.length ; i < count; i ++ ) {
    childIndex = childrenIndeces[i][0];
    goids.push(goNodes[childIndex]['i']);
  }

  return goids;
}*/


function getUniqueGoid(goids){
  var dic={},r=[];
  for(i=0;i<goids.length;i++){
    dic[goids[i][0]] = goids[i][1];
  }
  for(i in dic){
    arr = [i,dic[i]];
    r.push(arr);
  }
  return r;
}



//drawing the graph
function drawGraph(struct) {
  force
    .nodes(struct.nodes)
    .links(struct.links);

  var strokeColor = d3.scale.category10();
  var radiusScale = d3.scale.linear().domain([0,d3.max(struct.nodes,function(d) {return d.r;})]).range([5,10]);

  var link = go_chart.selectAll("line.link")
    .data(struct.links)
    .enter().append("line")
    .attr("class", "link")
    .style('stroke', function(d, i ) { return strokeColor(d.type);})
    .style("stroke-width", 2);
 
  var node = go_chart.selectAll("circle.node")
    .data(struct.nodes)
    .enter().append("circle")
    .attr("class", "node")
    .attr("r", function(d) { return radiusScale(d.r);})
    .style("fill", function(d,i) {
      if (d.is == "root") {
        return 'black';
      } else if (d.is == "target" || d.is == 'child') {
        return 'orange';
      } else {
        return 'lightblue';
      }
    })
    .call(force.drag)
    .on('click',function(d,i) {
      d.fixed=true;
      console.log("click");
    });
 
/*  var nodeText = go_chart.selectAll("text.node")
    .data(struct.nodes)
    .enter().append("text")
    .attr("class", "node")
    .attr("displayState",'y')
    //.style('display','none')
    .text(function(d) { return d.name;})
    .call(force.drag);*/



   
    var nodeText = go_chart.selectAll("text.node")
      .data(struct.nodes)
      .enter().append("text")
      .attr("class", "node")
      .text(function(d) { return d.name;})
      //.each(insertLinebreaks)
      .call(force.drag);

 
  force.on("tick", function(e) {
    link.attr("x1", function(d) { return d.source.x; })
      .attr("y1", function(d) { return d.source.y; })
      .attr("x2", function(d) { return d.target.x; })
      .attr("y2", function(d) { return d.target.y; });
 
    node.attr("cx", function(d) { return d.x; })
      .attr("cy", function(d) { return d.y; })

    nodeText.attr("x", function(d) { return d.x; })
      .attr("y", function(d) { return d.y + radiusScale(d.r) + 10; })
  });

  force.start();
}

var insertLinebreaks = function (d) {
          var el = d3.select(this);
          var words = d.name.split(' ');
          el.text('');
          var maxLength = 0;
          for (var i = 0; i < words.length; i++) {
              if(words[i].length > maxLength)
                  maxLength = words[i].length;
              if(words[i].length < 4)//if length is too short,connect to previous words to save space
              {
                words[i-1]=words[i-1]+" "+words[i];
                words[i]="";
              }
          }

          d.width = maxLength*9;
          d.height = words.length*9;

           for (var i = 0; i < words.length; i++) {
              if(words[i]!="")
              {
                var tspan = el.append('tspan').text(words[i]);
                tspan.attr('x', 0).attr('dy', '10')
                .attr("font-size", "3");
              }
          }

      };

//fixing and unfixing the nodes
function fixAllNodes() {
  for (var i = 0, count = force.nodes().length ; i < count ; i ++ ) {
    force.nodes()[i].fixed = true;
  }
}

function unfixAllNodes() {
  for (var i = 0, count = force.nodes().length ; i < count ; i ++ ) {
    if (force.nodes()[i]['is'] == "node") force.nodes()[i].fixed = false;
  }
}

//initializing statements


//recreate data structures and redraw graph with new nodes
function update(goid) {

  force = {};
  
  force = d3.layout.force()
  .charge(-3000)
  .linkDistance(50)
  .theta(0.2)
  .gravity(1)
  .size([width, height]);

  struct = parentGO(goid);
/*      immChildren = immediateChildren(goid);
  struct = struct.concat(immChildren);*/

  data = createD3Structure(struct,goid);

  go_chart.selectAll("line.link").remove();
  go_chart.selectAll("circle.node").remove();
  go_chart.selectAll("text.node").remove();

  drawGraph(data);
}

function getMaxPval(){
  var max = Number(go_inf[0].pVal);
  for(i=0;i<nodesSize;i++){
    if(Number(go_inf[i].pVal) > max){
      max = Number(go_inf[i].pVal);
    }
  }
  return max;
}

function sortGO_inf(){

  for(i=0;i<nodesSize;i++)
      pVal_sort[i] = Number(go_inf[i].pVal);

  function compare(a, b) {
  if (a < b) {
    return 1;
  }
  if (a > b) {
    return -1;
  }
  // a must be equal to b
  return 0;
  }


  pVal_sort.sort(compare);

}



function determineLabelSize(){
    var GOMinSize = 15;
    var GOMaxSize = 30;

    go_inf.map(function(d){
          d["size"] = getSize(d.pVal,GOMinSize,GOMaxSize);
    });

  function getSize(pVal,GOMinSize,GOMaxSize){
    for(i=0;i<seperated_points.length;i++){
        if((seperated_points[i] <= pVal) && (pVal <= seperated_points[i+1])){
                return (10-i)*(GOMaxSize-GOMinSize)/10+GOMinSize;
          }
      }

  }
}


function getColor(d,fill) {
/*    if(seperated_points[0] == d)
        return color_scheme[0];
*/
    /*console.log(d);*/
    for(i=0;i<seperated_points.length-1;i++){
      if((seperated_points[i] <= d)&&(d < seperated_points[i+1])){
          return color_scheme[i];
      }
    }
    if(seperated_points[seperated_points.length-1] == d)
      return color_scheme[9];


/*    if(seperated_points[seperated_points.length-1] == d)
        return color_scheme[9];*/
}

function createPvalLabelPanel(fill){
  titie = "p-value";
  panel.push(titie);
  panel.push("---------------------");
  //create panel
  for(i=0;i<seperated_points.length-1;i++)
  {
    from = seperated_points[i];
    to = seperated_points[i+1];
    panel.push('<i style="background:' + getColor(from,fill) + '" ></i> ' +
                from.toString().substring(0, 5) +" - "+ to.toString().substring(0, 5)
    );
  }

  panelLabel = panel.join('<br>');
  $('#pval-label').append(panelLabel);

}

function transformIndex(index)
{
  /*console.log(array_order);*/

  return (index < nodesSize)?
    array_order.indexOf(parseInt(index)):getIndexFormClusterNode(parseInt(index));

  function getIndexFormClusterNode(index){
    var nodeIndex = 0;
      clusterHierNodesStatus.map(function(d,i){
         if(d["index"] == index)
              nodeIndex = i;
      });
      return nodeIndex;
  }
    
}

function getMemKey(MenCache){
  var MemKeys = [];

  for (var i in MenCache){  
      MemKeys.push(i);
  }

  return MemKeys;
}

function transformPolarCoordinatesTORectangleCoordinates(radius,angel){
  x = radius*Math.sin(angel);
  y = radius*Math.cos(angel);
  return [x,y];
}

function darwClusterNodes(){


  clusterHierNodeView = clusterHierNodeView.data(clusterHierNodesStatus,function(d){return d.index;});

  clusterHierNodeView.enter().append("g");

  clusterHierNodeView.exit().remove();

  clusterHierNodeView.transition().duration(500).attrTween("transform", tween).attr("class","clusterNodeView");

  clusterHierNodeView.selectAll("circle").remove();

  circle = clusterHierNodeView.append("circle");

  circle.attr("r",function(d){
    return 6;
  }).style("fill", function(d){
    return (d["status"]=="Expanded")?"rgb(94, 141, 141)":"red";
  }).on("click",onClusterNodeClick)
  .append("title").text(function(d, i) {return "Click to cluster"});




  clusterHierNodeTextView.selectAll("text").remove();

  clusterHierNodeTextView = clusterHierNodeTextView.data(clusterHierNodesStatus,function(d){return d.index;});

  clusterHierNodeTextView.enter().append("g");

  clusterHierNodeTextView.exit().remove();

  clusterHierNodeTextView.transition().duration(500).attrTween("transform", tween).attr("class","clusterText");
  
  textBackground = clusterHierNodeTextView.append("g");

  text = textBackground.attr("transform", function(d) {
     return "rotate("  + (90 - d.angle * 180 / Math.PI) +")" + "translate(" + 0 + ",5)";
   }).append("text");

  text.attr('height', 'auto')
    .attr('text-anchor', 'middle')
    .text(function(d){
      return d["numOfOverlappingGenes"];
    }).attr('style', 'text-align:center;padding:2px;margin:2px;fill:white')
    .on("click",onClusterNodeClick)
    .append("title").text(function(d, i) {return "Click to cluster"});

  function tween(d, i, a) {
    var interpolate;
    var str;
    
    str = "rotate(" + (d["angle"] * 180 / Math.PI - 90) + ")"
     + "translate(" + ( d["radius"] + 5 ) + ",0)";
    interpolate = d3.interpolate(a,str);


    return function(t) {
        return interpolate(t);
      };
  }
}

function drawArc(clusterHierData){

  var clusterHierDataFiltered = [];
  var arcData = [];

  clusterHierData.map(function(d){
      if(d!=undefined)
        clusterHierDataFiltered.push(d);
  });
  
  clusterHierDataFiltered.map(function(d){
        var firstNodeIndex = d[0];
        var secondNodeIndex = d[1];
/*
        console.log(firstNodeIndex);
        console.log(secondNodeIndex);*/

        if(firstNodeIndex < nodesSize){
          firstNodeIndex = transformIndex(firstNodeIndex);
          
          firstNode = chordGroupsNodePosition[firstNodeIndex];
        }else{
          firstNodeIndex = transformIndex(firstNodeIndex);
          firstNode = clusterHierNodesStatus[firstNodeIndex];
        }
         if(secondNodeIndex < nodesSize){
          secondNodeIndex = transformIndex(secondNodeIndex);
          
          secondNode = chordGroupsNodePosition[secondNodeIndex];
        }else{
          secondNodeIndex = transformIndex(secondNodeIndex);
          secondNode = clusterHierNodesStatus[secondNodeIndex];
        }
  
    /*    console.log(firstNode);
        console.log(secondNode);*/
        if(firstNode.angle-secondNode.angle>Math.PI){
            startAngle = secondNode.angle;
            endAngle = firstNode.angle;
        }else{
          startAngle = firstNode.angle;
          endAngle = secondNode.angle;
        }

        if(firstNode.radius > secondNode.radius){
          radius = firstNode.radius;
        }
        else{
          radius = secondNode.radius;
        }

        arcData.push({"index":d[0],"radius":radius,"startAngle":startAngle,"endAngle":endAngle});

    });
        
  clusterArc = clusterArc.data(arcData,function(d){return d.index;});

  clusterArc.enter().append("path");

  clusterArc.exit().remove();

  clusterArc
    .style("fill", "green")
    .attr("class","clusterArc")
    .transition()
    .duration(500)
    .attr("d",function(d){
        
        var arc = d3.svg.arc()
          .innerRadius(d["radius"]+10)
          .outerRadius(d["radius"]+11)
          .startAngle(d["startAngle"])
          .endAngle(d["endAngle"]);

        return arc();
    });
}


 

function drawLine(clusterHierData){
  var clusterHierDataFiltered = [];
  var line1Data = [];
  var line2Data = [];
  
  clusterHierData.map(function(d){
      if(d!=undefined)
        clusterHierDataFiltered.push(d);
  })

  clusterHierDataFiltered.map(function(d){
        firstNodeIndex = d[0];
        secondNodeIndex = d[1];
        clusterNodeIndex = d[2];
        
        if(firstNodeIndex < nodesSize){
          firstNodeIndex = transformIndex(firstNodeIndex);
          firstNode = chordGroupsNodePosition[firstNodeIndex];
        }else{
          firstNodeIndex = transformIndex(firstNodeIndex);
          firstNode = clusterHierNodesStatus[firstNodeIndex];
        }

        if(secondNodeIndex < nodesSize){
          secondNodeIndex = transformIndex(secondNodeIndex);
          secondNode = chordGroupsNodePosition[secondNodeIndex];
        }else{
          secondNodeIndex = transformIndex(secondNodeIndex);
          secondNode = clusterHierNodesStatus[secondNodeIndex];
        }

        if(clusterNodeIndex < nodesSize){
          clusterNodeIndex = transformIndex(clusterNodeIndex);
          clusterNode = chordGroupsNodePosition[clusterNodeIndex];
        }else{
          clusterNodeIndex = transformIndex(clusterNodeIndex);
          clusterNode = clusterHierNodesStatus[clusterNodeIndex];
        }


        firstLineAngle = firstNode.angle;
        firstLineInnerPosition = firstNode.radius;
        firstLineOuterPosition = clusterNode.radius;

        secondLineAngle = secondNode.angle;
        secondLineInnerPosition = secondNode.radius;
        secondLineOuterPosition = clusterNode.radius;

        line1Data.push({"index":d[0],"LineInnerPosition":firstLineInnerPosition,"LineOuterPosition":firstLineOuterPosition,"LineAngle":firstLineAngle});

        line2Data.push({"index":d[1],"LineInnerPosition":secondLineInnerPosition,"LineOuterPosition":secondLineOuterPosition,"LineAngle":secondLineAngle});
  });
  
  clusterLine1 = clusterLine1.data(line1Data,function(d){return d.index;});

  clusterLine1.enter().append("path");

  clusterLine1.exit().remove();

  clusterLine1
    .style("fill", "green")
    .attr("class","clusterLine")
    .transition()
    .duration(500)
    .attr("d",function(d){
        
        var firstLine = d3.svg.arc()
        .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+5:d["LineInnerPosition"])
        .outerRadius(d["LineOuterPosition"]+5)
        .startAngle(d["LineAngle"])
        .endAngle(d["LineAngle"]+0.002);

        return firstLine();
    });

  clusterLine2 = clusterLine2.data(line2Data,function(d){return d.index;});

  clusterLine2.enter().append("path");

  clusterLine2.exit().remove();

  clusterLine2
    .style("fill", "green")
    .attr("class","clusterLine")
    .transition()
    .duration(500)
    .attr("d",function(d){

      var secondLine = d3.svg.arc()
        .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+5:d["LineInnerPosition"])
        .outerRadius(d["LineOuterPosition"]+5)
        .startAngle(d["LineAngle"])
        .endAngle(d["LineAngle"]+0.002);

        return secondLine();
    });
}

function drawHierarchicalClustering(clusterHierData){
    
    drawArc(clusterHierData);
    drawLine(clusterHierData);
    darwClusterNodes();
    
}


function redraw(transition) {

        (transition ? svg.transition() : svg)
          .attr("transform", "translate(" + zoom.translate() + ") scale(" + zoom.scale() + ")");

        if((zoom.scale()<0.8)&&(0.8<zoomLevel<0.9)){
              textBackground.attr("visibility","hidden");
          }

          if(zoom.scale()>=0.7)
          {
              textBackground.attr("visibility","visible");
          }

          zoomLevel = zoom.scale();
    
}


function getHierNodes(index){

  var leafNodes = [];
  var clusterNodesLevel = [];

  recursiveGetNodes(index-maxNodeIndex);

  function recursiveGetNodes(pos){

      var leftNode = clusterHierData[pos][0];
      var rightNode = clusterHierData[pos][1];

      if(leftNode < maxNodeIndex){
        leafNodes.push({"level":pos,"nodeId":leftNode});
      }else{
        recursiveGetNodes(leftNode-maxNodeIndex);
      }

      if(rightNode < maxNodeIndex){
        leafNodes.push({"level":pos,"nodeId":rightNode});
      }else{
        recursiveGetNodes(rightNode-maxNodeIndex);
      }

      clusterNodesLevel.push(pos);
  }

  return {"leafNodes":leafNodes,"clusterNodesLevel":clusterNodesLevel};
}

function getLeafNodesPosition(nodes){

    var nodePositions = [];
    nodes.map(function(d){
        nodePositions.push(transformIndex(d["nodeId"]));
    });

    return nodePositions;
}

function getClusterNodeLevel(nodes){
  var clusterNodesLevel = [];
  nodes.map(function(d){
    clusterNodesLevel.push(d["level"]);
  });

  return clusterNodesLevel.unique();
}

/*function getClusterNodesPosition(nodes){

    console.log(nodes[nodes.length-1]["level"]);
    var nodesPositions = [];

    
    return nodesPositions;
}*/

Array.prototype.unique = function() {
    var o = {}, i, l = this.length, r = [];
    for(i=0; i<l;i+=1) o[this[i]] = this[i];
    for(i in o) r.push(o[i]);
    return r;
};

function getClusterGenes(nodePositions){
  var clustergenes = [];

  nodePositions.map(function(d){ 
        go_inf[d].genes.map(function(d,i){
             clustergenes.push(d);
         });
  });
  
  return clustergenes.unique();
}


function calculateNewPvalue(removeGO){
  var pVal = 0;
  removeGO.map(function(d,i){

      pVal += parseFloat(d["pVal"]);
  });
  pVal = pVal/removeGO.length;
  return pVal;
}

function collapseNodeSet(nodePositions,genes,nodeBeingClicked){

  //console.log(go_inf);
  var removeGO = go_inf.splice(nodePositions[0],nodePositions.length);
  //console.log(removeGO);

 //console.log(go_inf);
  var GO_id = [];
  var GO_name = [];
  removeGO.map(function(d){
      GO_id.push(d["GO_id"]);
      GO_name.push(d["GO_name"]);
  });

  var pVal = calculateNewPvalue(removeGO);

  newGONode = {"GO_id":GO_id,"GO_name":GO_name,"count":genes.length,"pVal":pVal,"genes":genes,"nodeBeingClicked":nodeBeingClicked};

  go_inf.splice(nodePositions[0],0,newGONode);
  //console.log(go_inf);

  return removeGO;
}

function expandNodeSet(nodeBeingClicked,removeGO){
  var position;

  go_inf.map(function(d,i){
      if(d["nodeBeingClicked"] == nodeBeingClicked)
      {
        position = i;
      }
  });
  go_inf.splice(position,1);

  removeGO.map(function(d,i){
      go_inf.splice(position+i,0,d);
  });
  
}


function updateMatrix(){
  var newMatrix = [];
  var newArray = [];

  go_inf.map(function(d1,i){
      newArray = [];
      go_inf.map(function(d2,j){
        var numOfCommomGenes = 0;
        if(i != j)
          var numOfCommomGenes = getCommonGeneSize(d1["genes"],d2["genes"]);
        newArray.push(numOfCommomGenes);
      });
      newMatrix.push(newArray);
  });
  matrix = newMatrix;
  /*console.log(matrix);*/


  function getCommonGeneSize(genes1,genes2){
      var size = 0;
      genes1.map(function(d){
          if(genes2.indexOf(d)!=-1)
           size++;
      });
      return size;
  }
  return newMatrix;
}


function updateLayout(newMatrix){

    chord = chordMatrix.matrix(newMatrix);

    groupElements = groupElements.data(chord.groups());

    groupElements.enter().append("svg:path");

    groupElements.exit().remove();

    groupLayout =groupElements
     .style("fill", function(d) { return getColor(Number(go_inf[d.index].pVal),fill); })
     .style("stroke", function(d) { return getColor(Number(go_inf[d.index].pVal),fill); })
     .attr("d", d3.svg.arc().innerRadius(r0).outerRadius(r1)) //important
     .attr("id","chordGroups")
     .on("mouseover", mouseover_group);
/*     .append("title").text(function(d, i) {
      console.log(go_inf[d.index].GO_id);
      return go_inf[d.index].GO_id + "\t"+ go_inf[d.index].GO_name+ " , Num of genes:"+ go_inf[d.index].count + ", p value: " + go_inf[d.index].pVal ;
    });*/

    chordElements = chordElements.data(chord.chords());

    chordElements.enter().append("svg:path");

    chordElements.exit().remove();

    chordLayout = chordElements
      .attr("class", "chord")
      //.style("fill", function(d) { return fill(d.source.index%10); })
      //.style("opacity",".5")
      .attr("d", d3.svg.chord().radius(r0))
      .attr("id","chordChords")
      .on("mouseover", mouseover_chord);

    chordLayout.classed("fade", function(p) {
      return 1;
    });

}



function getMinValuePosition(nodePostions){
  var position = 0;
  var min = array_order[nodePostions[0]];

  nodePostions.map(function(d,i){
      if(array_order[d] < min){
        min = array_order[d];
        position = i;
      }
  });

  return position;

}


function calculateAClusterNodePosition(firstNodeIndex,secondNodeIndex,status,index,nodeBeingClicked,clusterNodesRadius,collapsedNodeRadius){
/*  console.log(array_order);
  console.log("nodesSize:" + nodesSize)
  console.log(firstNodeIndex);
  console.log(secondNodeIndex);*/

  //get the node object
  if(firstNodeIndex < nodesSize){
    firstNodeIndex = transformIndex(firstNodeIndex);
    firstNode = chordGroupsNodePosition[firstNodeIndex];
  }
  else{
    firstNodeIndex = transformIndex(firstNodeIndex);
    firstNode = clusterHierNodesStatus[firstNodeIndex];
  }

  if(secondNodeIndex < nodesSize){
    secondNodeIndex = transformIndex(secondNodeIndex);
    secondNode = chordGroupsNodePosition[secondNodeIndex];
  }
  else{
    secondNodeIndex = transformIndex(secondNodeIndex);
    secondNode = clusterHierNodesStatus[secondNodeIndex];
  }

/*  console.log("t:"+firstNodeIndex);
  console.log("t:"+secondNodeIndex);
  console.log(firstNode);
  console.log(secondNode);
*/
  //console.log(clusterNodesRadius);

  //calculate the position of new node(radius,angle)
  angle = (firstNode.angle+secondNode.angle)/2;

  if((nodeBeingClicked == index)&&(status=="collapse")){
    radius = collapsedNodeRadius;
  }else if((clusterNodesRadius!=undefined)&&(clusterNodesRadius.map(function(d){return d["nodeId"]}).indexOf(index)!=-1)){
    clusterNodesRadius.map(function(d){
      if(d["nodeId"]==index){
          radius = d["radius"];
      }
    });
  }else{
     radius = (firstNode.radius > secondNode.radius)?firstNode.radius+5:secondNode.radius+5;
  }
  return [angle,radius];
}

function createClusterNodesStatus(clusterHierData,nodeBeingClicked,status,collpasedNodes,clusterNodesRadius,collapsedNodeRadius){

/*  console.log(clusterHierData);*/
  for(i=0;i<clusterHierData.length;i++){
    if(clusterHierData[i]!=undefined){
        firstNode = clusterHierData[i][0];
        secondNode = clusterHierData[i][1];
        index = clusterHierData[i][2];
        numOfOverlappingGenes = clusterHierData[i][3];

        position = calculateAClusterNodePosition(firstNode,secondNode,status,index,nodeBeingClicked,clusterNodesRadius,collapsedNodeRadius);
        //create new cluster node based on the position
        if(nodeBeingClicked == index){
          if(status == "collapse"){
            nodeObj = {"index":index,"angle":position[0],"radius":position[1],"numOfOverlappingGenes":numOfOverlappingGenes,"status":"Collapsed"};
          }else{
            nodeObj = {"index":index,"angle":position[0],"radius":position[1],"numOfOverlappingGenes":numOfOverlappingGenes,"status":"Expanded"};
          }
        }
        else{
          if(collpasedNodes.indexOf(index.toString())!=-1){
            nodeObj = {"index":index,"angle":position[0],"radius":position[1],"numOfOverlappingGenes":numOfOverlappingGenes,"status":"Collapsed"};
          }else{
            nodeObj = {"index":index,"angle":position[0],"radius":position[1],"numOfOverlappingGenes":numOfOverlappingGenes,"status":"Expanded"};
          }
        }
        
        clusterHierNodesStatus.push(nodeObj);
      }
    }
}


/*
*get rid of the merge node
*/
function removeNodesInArrayOrder(nodePositions){
  //remove nodes in array
      var newArrayOrder = [];
      var removeNodeInArray = [];
/*
      console.log("array_order = " + array_order);

      console.log("nodePositions:"+nodePositions);*/

      array_order.map(function(d,i){
          if(nodePositions.indexOf(i) == -1 ){
                  newArrayOrder.push(d);
          }
          else{
            if(nodePositions.indexOf(i) == nodePositions.length -1){
              newArrayOrder.push(d);
            }
            removeNodeInArray.push(d);
          }
       
      });
      old_array_order = array_order;
      array_order = newArrayOrder;
      
      /*console.log("array_order = " + array_order);*/

      return removeNodeInArray;
}

function addNodesToArrayOrder(removedNodes){
  var position;

  array_order.map(function(d,i){
    if(removedNodes.indexOf(d)!=-1){
      position = i;
    }
  });

  array_order.splice(position,1);

  removedNodes.map(function(d,i){
      array_order.splice(position+i,0,d);
  });

  /*console.log("recover array order:"+array_order);*/
}

function removeNodesInClusterData(nodePositions,clusterNodesLevel){
    /*console.log(clusterNodesLevel);*/
    var newNode = old_array_order[nodePositions[nodePositions.length-1]];
    /*console.log("newNode:"+newNode);*/
    var newClusterHierData = [];

    clusterHierData.map(function(d,i){
        if(clusterNodesLevel.indexOf(i) == -1){
            newClusterHierData.push(d);
        }else{
            if(clusterNodesLevel.indexOf(i) == clusterNodesLevel.length-1)
            {
                newClusterHierData.push([newNode,newNode,d[2]]);
            }else{
                newClusterHierData.push(undefined);
            }
        }
    });    

    /*console.log(newClusterHierData);*/

    return newClusterHierData;
}

function addNodesInClusterData(removedClusterHierData){

/*  clusterHierData.map(function(d,i){
      if(clusterNodesLevel.indexOf(i)!=-1){
          clusterHierData[i] = removedClusterHierData
      }
  });*/

  removedClusterHierData.map(function(d,i){
      clusterHierData[d[2]-maxNodeIndex] = d;
  });

  return clusterHierData;
}

function getRemoveClusterData(clusterNodesLevel){
  var removeHierData = [];

  clusterNodesLevel.map(function(d,i){
      removeHierData.push(clusterHierData[d]);
  });

  return removeHierData;
}

function collapseHierClustering(nodePositions,clusterNodesLevel,nodeBeingClicked,collpasedNodes,clusterNodesRadius){
    var collapsedNodeRadius;
    clusterHierNodesStatus.map(function(d){
      if(d["index"]==nodeBeingClicked)
        collapsedNodeRadius = d["radius"];
    });
     
    clusterHierNodesStatus = [];//empty clusterHierNodesStatus
    chordGroupsNodePosition = [];

    /*console.log(clusterNodesLevel);*/

    chord.groups().forEach(function(d,i){//recalculate angel,radius
      nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
      chordGroupsNodePosition.push(nodeObj);
    });

    var removeHierData = getRemoveClusterData(clusterNodesLevel);
    /*console.log(removeHierData);*/

    clusterHierData = removeNodesInClusterData(nodePositions,clusterNodesLevel);

    createClusterNodesStatus(clusterHierData,nodeBeingClicked,"collapse",collpasedNodes,clusterNodesRadius,collapsedNodeRadius);
    
    /*console.log("newclusternode:"+clusterHierData);*/
    if(labelOff==1)
      drawHierarchicalClustering(clusterHierData);

    return removeHierData;
}

function expandHierClustering(removedClusterHierData,clusterNodesLevel,nodeBeingClicked,collpasedNodes,clusterNodesRadius){
    clusterHierNodesStatus = [];//empty clusterHierNodesStatus
    chordGroupsNodePosition = [];

    chord.groups().forEach(function(d,i){//recalculate angel,radius
      nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
      chordGroupsNodePosition.push(nodeObj);
    });

    clusterHierData = addNodesInClusterData(removedClusterHierData,clusterNodesLevel);

    createClusterNodesStatus(clusterHierData,nodeBeingClicked,"expand",collpasedNodes,clusterNodesRadius);

    if(labelOff==1)
      drawHierarchicalClustering(clusterHierData);
}

function onClusterNodeClick(d,i){

/*  console.log("node being clicked");
  console.log(d);*/

  //transform the index to clusterlevel

  var nodeBeingClicked = d["index"];

  if(d["status"]=="Expanded"){
    /*console.log("node going to Collapse");*/
    var nodes = getHierNodes(nodeBeingClicked);
/*    console.log(nodes);*/

    var nodePositions = getLeafNodesPosition(nodes["leafNodes"]).unique();
/*    console.log("array_order:" + array_order);
    console.log(nodePositions);*/

    var clusterNodesLevel = nodes["clusterNodesLevel"];

    var genes = getClusterGenes(nodePositions);

    var removeGO = collapseNodeSet(nodePositions,genes,nodeBeingClicked);

    createInteractionGenes(go_inf);
    createGeneInfo();

    matrix = updateMatrix();
    updateLayout(matrix);

    var removeNodeInArray = removeNodesInArrayOrder(nodePositions);
    /*console.log("array_order:" + array_order);*/
    /*console.log("clusterNodeLevel:" + nodes["clusterNodeLevel"]);*/

    var removeHierData = collapseHierClustering(nodePositions,clusterNodesLevel,nodeBeingClicked,Object.keys(memCache),memCache["clusterNodesRadius"]);

    if(memCache["clusterNodesRadius"]!=undefined)
        var clusterNodesRadius = memCache["clusterNodesRadius"];
    else
        clusterNodesRadius = [];

    clusterHierNodesStatus.map(function(d){
        if(d["index"]==nodeBeingClicked){
          //console.log(d);
          clusterNodesRadius.push({"nodeId":d["index"],"radius":d["radius"]});
        }
    });

    memCache["clusterNodesRadius"] = clusterNodesRadius;

    //console.log(clusterNodesRadius);
    /*console.log(removeHierData);*/

    var nodeBeingMemorized = {"go_inf":removeGO,"array_order":removeNodeInArray,"clusterHierData":removeHierData,"clusterNodesLevel":clusterNodesLevel};

    memCache[nodeBeingClicked] = nodeBeingMemorized;

  }else{
    /*console.log("node going to expand");*/

    nodeBeingMemorized = memCache[nodeBeingClicked];

    expandNodeSet(nodeBeingClicked,nodeBeingMemorized["go_inf"]);

    createInteractionGenes(go_inf);

    createGeneInfo();

    matrix = updateMatrix();

    updateLayout(matrix);

    addNodesToArrayOrder(nodeBeingMemorized["array_order"]);

    var collapsedNodes = getMemKey(memCache);

    var clusterNodesLevel = nodeBeingMemorized["clusterNodesLevel"];

    var clusterNodesRadius = memCache["clusterNodesRadius"];

    var newclusterNodesRadius = [];

    clusterNodesRadius.map(function(d,i){
      if(d["nodeId"]!=nodeBeingClicked)
        newclusterNodesRadius.push(d);
    })

    clusterNodesRadius = newclusterNodesRadius;

    memCache["clusterNodesRadius"] = clusterNodesRadius;

    expandHierClustering(nodeBeingMemorized["clusterHierData"],clusterNodesLevel,nodeBeingClicked,collapsedNodes,clusterNodesRadius);

    delete memCache[nodeBeingClicked];
  }
}

//draw lines,arcs,node



function moveOutHierCluster(){
  var clusterNode = d3.selectAll(".clusterNode");

  clusterNode.transition().duration(500).attr("transform",function(d) {
     return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")"
         + "translate(" + (d.radius+10) + ",0)";
   });

  var clusterArc = d3.selectAll(".clusterArc");

  clusterArc.transition().duration(500).attr("d",function(d){
        var arc = d3.svg.arc()
          .innerRadius(d["radius"]+10)
          .outerRadius(d["radius"]+11)
          .startAngle(d["startAngle"])
          .endAngle(d["endAngle"]);
        return arc();
  });

  var clusterLine = d3.selectAll(".clusterLine");

  clusterLine.transition().duration(500).attr("d",function(d){
        var Line = d3.svg.arc()
        .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+5:d["LineInnerPosition"])
        .outerRadius(d["LineOuterPosition"]+5)
        .startAngle(d["LineAngle"])
        .endAngle(d["LineAngle"]+0.002);

        return Line();
  });

}

function moveInHierCluster(){
  var clusterNode = d3.selectAll(".clusterNode");

  clusterNode.transition().duration(500).attr("transform",function(d) {
     return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")"
         + "translate(" + (d.radius) + ",0)";
   });

  var clusterArc = d3.selectAll(".clusterArc");

  clusterArc.transition().duration(500).attr("d",function(d){
        var arc = d3.svg.arc()
          .innerRadius(d["radius"]+5)
          .outerRadius(d["radius"]+6)
          .startAngle(d["startAngle"])
          .endAngle(d["endAngle"]);
        return arc();
  });

  var clusterLine = d3.selectAll(".clusterLine");

  clusterLine.transition().duration(500).attr("d",function(d){
        var Line = d3.svg.arc()
        .innerRadius(d["LineInnerPosition"])
        .outerRadius(d["LineOuterPosition"])
        .startAngle(d["LineAngle"])
        .endAngle(d["LineAngle"]+0.002);

        return Line();
  });
}


function getLevelFromNumOfOverlappingGenes(numOfOverlappingGenes){
  var level = -1;
  clusterHierDataStatic.map(function(d,i){
      if(d[1]>=numOfOverlappingGenes)
        level = i;
  });

  return level;
}



function getClusterNodesIndexBeingSelected(level){

  var clusterDataLevel = [];
  if(level>=level_g){//collapse
    
    for(var i=level;i >= 0; i--){
      if(clusterHierData[i]!=undefined&&clusterDataLevel.indexOf(i)==-1){
        if(memCache[clusterHierData[i][2]]==undefined){
          var clusterNodesLevel = getHierNodes(clusterHierData[i][2])["clusterNodesLevel"];
          //console.log(clusterNodesLevel);
          clusterNodesLevel.map(function(d){return clusterDataLevel.push(d);});
          clusterDataLevel = clusterDataLevel.unique();
          clusterNodesIndex.push(clusterHierData[clusterNodesLevel[clusterNodesLevel.length-1]][2]);
        }
      }
      clusterNodesIndex = clusterNodesIndex.unique();


    }
    
  }else{//expand
    var index = level + nodesSize;
    var pos = [];
    clusterNodesIndex.reverse().map(function(d1,i){
      if(d1 >= index){
/*        clusterHierNodesStatus.map(function(d){
          if(d["index"] == d1){
            console.log(d1);
            onClusterNodeClick(d);
          }
        });*/
        pos.push(i);

        var clusterNode = getClusterNode(d1);
        onClusterNodeClick(clusterNode);
      }
    });

    clusterNodesIndex.splice(pos[0],pos.length);

    for(var i=level;i >= 0; i--){
      if(clusterHierData[i]!=undefined&&clusterDataLevel.indexOf(i)==-1){

        var clusterNodesLevel = getHierNodes(clusterHierData[i][2])["clusterNodesLevel"];
        //console.log(clusterNodesLevel);
        clusterNodesLevel.map(function(d){return clusterDataLevel.push(d);});
        clusterDataLevel = clusterDataLevel.unique();
        clusterNodesIndex.push(clusterHierData[clusterNodesLevel[clusterNodesLevel.length-1]][2]);

      }
    }
    clusterNodesIndex = clusterNodesIndex.unique();

  }

  if(clusterNodesIndex.length!=0){
    clusterNodesIndex.map(function(index,i){
        clusterHierNodesStatus.map(function(d){
          if(index == d["index"]){
            if(memCache[index]==undefined){
              onClusterNodeClick(d);
            }
          }
        });
        
    });
  }
  level_g = level;
}

function getClusterNode(index){
  var node;
  clusterHierNodesStatus.map(function(d){
    if(d["index"] == index)
      node = d;
  });
  return node;
}

function createGoHier(goid){
  
  go_chart = d3.select("#go_chart").append("svg")
    .attr("width", width)
    .attr("height", height);
  go_chart.append('rect')
    .style('fill','white')
    .style('stroke','gray')
    .attr('width',width)
    .attr('height',height)
    .attr('x',0)
    .attr('y',0);

  update(goid);
}


 function mouseover_group(d, i) {
    $('#content').empty();
    chordLayout.classed("fade", function(p) {
      return p.source.index != i
          && p.target.index != i;
    });
    var genes = "";

    go_inf[d.index].genes.forEach(function(d,i){

       var tmp = "<p>" + (i+1) +"."+ d + "</p>";
       genes+=tmp;
    });

    var str = "<p> GO_id: " + go_inf[d.index].GO_id + "<p> GO_Name: "+ go_inf[d.index].GO_name+ " <p> Num of genes: "+ go_inf[d.index].count + "<p> P-value: " + go_inf[d.index].pVal;


    var words = [];

/*    if(typeof go_inf[d.index].GO_name != "string"){
      go_inf[d.index].GO_name.map(function(d) {
          words.push({text: d, size: goLabelSize[d]});
      });

      $('#content').append("GO name:");
      d3.layout.cloud().size([450, 450])
        .words(words)
        .rotate(function() { return ~~(Math.random() * 5) * 10; })
        .font("Impact")
        .fontSize(function(d) { return d.size; })
        .on("end", draw)
        .start();

      function draw(words) {
        var fill = d3.scale.category20();

        d3.select("#content").append("svg")
            .attr("width", 450)
            .attr("height", 450)
          .append("g")
            .attr("transform", "translate(225,225)")
          .selectAll("text")
            .data(words)
          .enter().append("text")
            .style("font-size", function(d) { return d.size + "px"; })
            .style("font-family", "Impact")
            .style("fill", function(d, i) { return fill(i); })
            .attr("text-anchor", "middle")
            .attr("transform", function(d) {
              return "translate(" + [d.x, d.y] + ")rotate(" + d.rotate + ")";
            })
            .text(function(d) {console.log(d.text); return d.text; });
      }

      $('#content').append("<p>Genes:\n"+genes);
    }else{*/
      $('#content').append(str + "<p>Genes:\n"+genes);
    //}

    if(typeof go_inf[d.index].GO_id == "string"){
      $('#content').append("<div id=\"go_chart\"></div>");
      createGoHier(go_inf[d.index].GO_id);
    }
}

 function mouseover_chord(d, i) {
    $('#content').empty();
    var index = d.source.index+"-"+d.source.subindex;

    var str = "";

    dic[index].split("|").forEach(function(d,i){
       var tmp = "<p>" + (i+1) +"."+ d + "</p>";
       str+=tmp;
    });

    $('#content').append("<p>Overlapping genes between " + go_inf[d.source.index].GO_id +" and " +  go_inf[d.source.subindex].GO_id+":\n</p>"+str);

}




function createOnClick(num_array){
  while(num_array.length!=0){
    var num = num_array.pop();
    var term = "#GO_button_"+num;
    $content = $('#content').on('click',"#GO_button_"+num, function(d) {
            var num = d.toElement.id.split("_")[2];
            chordLayout.classed("fade", function(p) {
              return p.source.index != num
                  && p.target.index != num;
            });
    });
  }
    
    $(".dropbtn").click(function () {

      $header = $(this);
      //getting the next element
      $content = $header.next();

      //open up the content needed - toggle the slide- if visible, slide up, if not slidedown.
      $content.slideToggle(500, function () {
          //execute this after slideToggle is done
          //change text of header based on visibility of content div
  /*        $header.text(function () {
              //change text based on condition
              return $content.is(":visible") ? "Collapse" : "Expand";
          });*/
      });
    });

}

function getGenesFromACluster(num){
  var genes = "";
  go_inf[num].genes.forEach(function(d,i){

       var tmp = "<p>" + (i+1) +"."+ d + "</p>";
       genes+=tmp;
    });
  return genes;
}

function renderGOTerm(go_num_array){
  var GOterms="";
  for(i=0;i<go_num_array.length;i++){
    num = go_num_array[i];
    Go = go_inf[num].GO_id + " " + go_inf[num].GO_name;
    var line = "<button class=\"dropbtn\""+ "id=GO_button_"+ num +">" + Go + "</button>";
    line += "<div class=\"Go_content\">" + "<p> P-value: " + go_inf[num].pVal + " <p> Num of genes: "+ go_inf[num].count +
    "<p>Genes: " +  getGenesFromACluster(num) + "</div>";
    GOterms += line; 
   }
  return GOterms;
}

function getGoIndex(go_id){
  for(i in go_inf){
      if(go_inf[i].GO_id == go_id)
        return i;
  }
}

function isNum(str){
  if(str=="") return false;
  arr = str.split('');
  for(i=0;i<arr.length;i++)
  {
    if(arr[i]<'0'||arr[i]>'9')
      return false;
  }
  return true;
}

function sortGOcontent(go_contents){

  function compare(a, b) {
    pa = Number(go_inf[a].pVal);
    pb = Number(go_inf[b].pVal);
    if (pa > pb) {
      return 1;
    }
    if (pa < pb) {
      return -1;
    }
    // a must be equal to b
    return 0;
  }
  go_contents.sort(compare);
}



function resetVis(){
    //text.transition().duration(500).text(function(d){return go_inf[d.index].GO_id;}).attr("x",8);
    groupLayout.transition().duration(500).attr("d",d3.svg.arc().innerRadius(r0).outerRadius(r1));
    chordLayout.transition().duration(500).attr("d",d3.svg.chord().radius(r0));
    //moveInHierCluster();
}

function refreshDetailPanel(){
  $('#content').empty();

  function isGO(searchTerm){
    if(isNum(searchTerm)||searchTerm.split(":")[0].toUpperCase()=="GO")
      return true;

    for(var i=0;i<go_inf.length;i++){
      if(go_inf[i]["GO_name"] == searchTerm.toLowerCase()){
        return true;
      }
    }


    return false;
  }
  
  var searchTerm = $('#filter').val(); 

  if(isGO(searchTerm)){//for GO
    resetVis();
    var i;

    if(isNum(searchTerm)){
      searchTerm = "GO:"+searchTerm;
      i = getGoIndex(searchTerm.toUpperCase());

    }else if(searchTerm.split(":")[0].toUpperCase()=="GO"){

      i = getGoIndex(searchTerm.toUpperCase());
    }else{

      for(var j=0;j<go_inf.length;j++){
        if(go_inf[j]["GO_name"] == searchTerm.toLowerCase()){
          i = j;
        }
      }
    }
    

    if(i != undefined){
    
      str = "<p> GO_id: " + go_inf[i].GO_id + "<p> GO_Name: "+ go_inf[i].GO_name+ " <p> Num of genes: "+ go_inf[i].count + "<p> P-value: " + go_inf[i].pVal ;
      genes = getGenesFromACluster(i);
      $('#content').append(str + "<p>Genes:\n"+genes+"<div id=\"go_chart\"></div>");
      createGoHier(go_inf[i].GO_id);

      chordLayout.classed("fade", function(p) {
              return p.source.index != i
                  && p.target.index != i;
            });

      groupLayout.transition().attr("d",
           d3.svg.arc().innerRadius(function(d){return (d.index!=i)? r0:r0+5;}).outerRadius(
            function(d){return (d.index!=i)? r1:r1+5;}));

      chordLayout.transition().attr("d",d3.svg.chord().radius(function(d){return (d.index!=i)? r0:r0+5;}));

      
      //moveHierCluster();
    }
  }else{//for genes
  var indexs = gene_info[searchTerm];

  
  if(indexs!=undefined)
  {

      var index_array = indexs.split(";");
      var go_contents = [];
      var num_array= [];

      for(i=0;i<index_array.length;i++){

        var num = parseInt(index_array[i]);
        num_array.push(num);
        
        go_contents.push(num);
      }

      sortGOcontent(go_contents);

      $('#content').append(renderGOTerm(go_contents));

      createOnClick(num_array);

      var arrayList=[];
      for(i in index_array){
        var num = parseInt(index_array[i]);
        arrayList.push(num);
      }

 /*     text.transition().duration(500).text(function(d){
        if(arrayList.indexOf(d.index)==-1)
          return go_inf[d.index].GO_id;
        else
          return "["+go_inf[d.index].GO_id+"]";
        
        
      }).attr("x",function(d){
        if(arrayList.indexOf(d.index)!=-1){
            return d.angle > Math.PI ? 1:16;
          }else
          {
            return 8;
          }
         });*/

      groupLayout.transition().attr("d",
           d3.svg.arc().innerRadius(function(d){return (arrayList.indexOf(d.index)==-1)? r0:r0+5;}).outerRadius(
            function(d,i){return (arrayList.indexOf(d.index)==-1)? r1:r1+5;}));

      chordLayout.transition().attr("d",d3.svg.chord().radius(function(d){return (arrayList.indexOf(d.index)==-1)? r0:r0+5;}));

   /*   if(arrayList.length!=0)
        moveOutHierCluster();*/
  }
  else{
    resetVis();
   
  }
}
}


function drawLable(){

  labelOff = 0;

  clusterArc.style("display","none");
  clusterLine1.style("display","none");
  clusterLine2.style("display","none");
  clusterHierNodeView.style("display","none");
  clusterHierNodeTextView.style("display","none");
  groupText.style("display","")

  updateLabel();
/*      .text(function(d) {return d.text});*/
/*   createSimpleClouds();

   function createSimpleClouds(words,size,$text){
      words.map(function(d){
        d3.select($text)
          .append("svg:text")
           .attr("x", 8)
           .attr("dy", ".45em")
           .attr("text-anchor", function(d) {
             return d.angle > Math.PI ? "end" : null;
           })
           .attr("transform", function(d) {
             return d.angle > Math.PI ? "rotate(180)translate(-16)" : null;
           })
           .text(d.text);
      });
   }*/
/*   var angle = 0;
    var textPos = [];
   var sizes = [];

   for(var i =0;i<textPos.length;i++){
    if(i == textPos.length-1){
      angleWidth = 2*Math.PI - textPos[i].angle;
      console.log(angleWidth);
    }else{
      angleWidth = textPos[i+1].angle - textPos[i].angle;
    }

      var width = 2*angleWidth*r1;
      var height = 2*angleWidth*r1;
      sizes.push([width,height]);
   }

   text.each(function(d){
      var gos = go_inf[d.index].GO_id;
      var words = [];
      if(typeof gos == "string"){
          words.push({"text":gos,"size":10});
      }else{
        
         gos.map(function(d){
            words.push({"text":d,"size":10});
        });
      }

      createClouds(words,sizes[d.index],this);
   })*/


 /* if(clouds.length>0){

      cloudsText = svg.append("svg:g")
       .selectAll("g")
         .data(clouds)
        .enter().append("svg:g")
        .attr("transform", function(d) {
            return "rotate(" + (d.startAngle * 180 / Math.PI - 90) + ")"
         + "translate(" + r1 + ",0)";
        });

        cloudsText.each(function(d){
          var width = (d.endAngle - d.startAngle)*r1;
          createClouds(d.words,width,this);
        });
        
  } 
   */

/*   function createClouds(word,size,$text){
      var words=[];
      d3.layout.cloud().size([size*2,size])
        .words(word)
        .rotate(function() { return ~~(Math.random() * 2) * 90; })
        .font("Impact")
        .fontSize(function(d) {words.push(d); return d.size;})
        .on("end", function() {draw(words,$text,size)})
        .start();
   }

   
   function draw(words,$text,size) {

    console.log(words);

      d3.select($text)
      .append("svg")
        .attr("width", size*2)
        .attr("height", size)
        .append("g")
        .attr("transform", function(d){

          return "translate("+(size)+","+(size/2)+")"
        })
      .selectAll("text")
        .data(words)
      .enter().append("text")
        .style("font-size", function(d) { return d.size + "px"; })
        .style("font-family", "Impact")
        .style("fill", function(d, i) { return fill(i); })
         .attr("text-anchor", function(d) {
           return d.startAngle > Math.PI ? "end" : null;
         })
        .attr("transform", function(d) {
          return  "translate(" + [d.x/2, d.y/2] + ")rotate(" + d.rotate+180 + ")";
        })
        .text(function(d) { return d.text; });
    }*/
}


function updateLabel(){
  groupText.remove();

  groupText = svg.append("svg:g")
        .selectAll("g")
        .data(chord.groups)
        .enter().append("svg:g")
        .attr("transform", function(d) {
        return "rotate(" + ((d.startAngle+d.endAngle)/2 * 180 / Math.PI - 90) + ")"
           + "translate(" + r1 + ",0)";
        }).append("svg:text")
        .attr("x", 8)
        .attr("dy", ".45em")
        .attr("text-anchor", function(d) {
         return (d.startAngle+d.endAngle)/2 > Math.PI ? "end" : null;
        })
        .attr("transform", function(d) {
         return (d.startAngle+d.endAngle)/2 > Math.PI ? "rotate(180)translate(-16)" : null;
        })
        .text(function(d) {return (typeof go_inf[d.index].GO_id == "string")? go_inf[d.index].GO_name:getMaxLabel(go_inf[d.index].GO_name)+"+"});



/*  var groups = [];
  var clouds = [];
  groupText.data().map(function(d,i){
   var words = [];

    if(typeof go_inf[i].GO_id == "string"){
        groups.push({"endAngle":d.endAngle,"startAngle":d.startAngle,"value":d.value,"text":go_inf[i].GO_id});
    }else{
      go_inf[i].GO_id.map(function(data){
        var str = data.toString().split(",");
        str.map(function(d){
             words.push({"text":d,"size":Math.ceil(goLabelSize[d])});
        });
      });
      clouds.push({"endAngle":d.endAngle,"startAngle":d.startAngle,"value":d.value,"words":words});
    }

   
  });*/
}

function getMaxLabel(d){
  var position = 0;
  var maxSize = 0;
  d.map(function(d,i){
    if(maxSize < goLabelSize[d]){
      maxSize = goLabelSize[d];
      position = i;
    }
  });

  return d[position];
}


/*  var textClouds = [];

  groups.map(function(d,i){
    
    var width = (d.endAngle - d.startAngle)*r1;
    var size = 0;
    d.words.map(function(d){
      size += d.size;
    });

    k = size/width;

    var size2 = 0;
    d.words.map(function(word,j){
        if(j==0){
          var angle = d.startAngle;
        }else{
          var angle = k*size2/r1 + d.startAngle;
        }
        textClouds.push({"text":word.text,"size":word.size,"angle":angle});
        size2 += word.size;
    });


  });*/

function drawHierCluster(){
  labelOff = 1;
  clusterArc.style("display","");
  clusterLine1.style("display","");
  clusterLine2.style("display","");
  clusterHierNodeView.style("display","");
  clusterHierNodeTextView.style("display","");


  groupText.style("display","none");
  //cloudsText.style("display","none");

  drawHierarchicalClustering(clusterHierData);


  
}

var force;

var nodesSize = size;
var maxNodeIndex = size;


var old_array_order = [];
var clusterHierDataFiltered = [];


var clusterHierDataStatic = [];
clusterHierData.map(function(d){
  clusterHierDataStatic.push([d[2],d[3]]);
});

var maxNumOfOverlappingGens = clusterHierData[0][3];
var minNumOfOverlappingGens = clusterHierData[clusterHierData.length-1][3];

var memCache = {};
var gene_info = {};

var labelOff=0;

var clusterHierNodesStatus = [];//for storing nodes generated Hierarchical clustering
var chordGroupsNodePosition = [];


go_inf.forEach(function(d,i){
  go_inf[i].genes = d.genes.split("|");
});

createGeneInfo();


var dic = {};// for storing gene intersection information
createInteractionGenes(go_inf);


/*GOFilter();
updateMatrix();*/

var chordMatrix = d3.layout.chord()
   .padding(.03)
   .sortSubgroups(d3.descending);
   
var chord = chordMatrix.matrix(matrix);
 
var w = $(window).width(),
     h = $(window).height(),
     r0 = Math.min(w, h) * .25,
     r1 = r0 * 1.1;

var width = 400;
var height = 600;
var force = d3.layout.force()
  .charge(-3000)
  .linkDistance(50)
  .theta(0.2)
  .gravity(1)
  .size([width, height]);
var svg ;
var go_chart;


var maxpVal = getMaxPval();

var minpVal = getMinPval();

var pVal_sort=[];

sortGO_inf();

var range = maxpVal-minpVal;
var step = range/10;
var seperated_points = [];
var panel=[];


for(i=0;i<11;i++){
  seperated_points.push(minpVal+step*i);
}

var fill = d3.scale.ordinal()
     .domain(d3.range(12).reverse())
     .rangeRoundPoints([1,255]);

var color_scheme = ["#1249C9","#0F399B","#0A2B76","#0F2147","#2C3645","#82733D","#E1C03B","#E6AE29","#F2AB1C","#FFAB00"].reverse();


createPvalLabelPanel(fill);

var zoom = d3.behavior.zoom().translate([w/2, h/2]);


var svg = d3.select("#chart")
   .append("svg:svg")
     .attr("class","main_vis")
     .attr("width", w)
     .attr("height", h)
     .call(zoom.on("zoom", redraw))
   .append("svg:g")
     .attr("id", "circle")
     .attr("transform", "translate(" + w / 2 + "," + h / 2 + ")");


svg.append("circle")
  .attr("r", r1);


 //store the nodes for hierarchical clustering visualzaition
 chord.groups().forEach(function(d,i){
  nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
  chordGroupsNodePosition.push(nodeObj);
 });

 createClusterNodesStatus(clusterHierData,[],"",[]);

var clusterArc = svg.append("svg:g")
 .selectAll("g")
   .data(clusterHierData)
 .enter().append("path");

var clusterLine1 = svg.append("svg:g")
 .selectAll("g")
   .data(clusterHierData)
 .enter().append("path");

var clusterLine2 = svg.append("svg:g")
 .selectAll("g")
   .data(clusterHierData)
 .enter().append("path");


var clusterHierNodeView = svg.append("svg:g")
 .selectAll("g")
   .data(clusterHierNodesStatus)
 .enter().append("svg:g");

var circle = clusterHierNodeView.attr("transform", function(d) {
     return "rotate(" + (d.angle * 180 / Math.PI - 90) + ")"
         + "translate(" + (d.radius+5) + ",0)";
   })
 .attr("id","clusterNode")
 .append("circle")
 .append("title").text(function(d, i) {return "Click to cluster"});

var clusterHierNodeTextView = svg.append("svg:g")
 .selectAll("g")
  .data(clusterHierNodesStatus)
 .enter().append("svg:g");

drawHierarchicalClustering(clusterHierData);



var groupElements = svg.append("svg:g")
 .selectAll("path")
   .data(chord.groups)
   .enter().append("svg:path");

var groupLayout = groupElements
   .style("fill", function(d) { return getColor(Number(go_inf[d.index].pVal),fill); })
   .style("stroke", function(d) { return getColor(Number(go_inf[d.index].pVal),fill); })
   .attr("d", d3.svg.arc().innerRadius(r0).outerRadius(r1)) //important
   .attr("id","chordGroups")
   .on("mouseover", mouseover_group);


var goLabelSize = {};

determineLabelSize();

go_inf.map(function(d){
  if(goLabelSize[d.GO_name] == undefined){
    goLabelSize[d.GO_name] = d.size;
  }
})



var chordElements = svg.append("svg:g")
      .selectAll("path")
      .data(chord.chords)
      .enter().append("svg:path");

var chordLayout = chordElements
      .attr("class", "chord")
      //.style("fill", function(d) { return fill(d.source.index%10); })
      //.style("opacity",".5")
      .attr("d", d3.svg.chord().radius(r0))
      .attr("id","chordChords")
      .on("mouseover", mouseover_chord);

/*chordLayout.append("title").text(function(d, i) {

    var index = d.source.index+"-"+d.source.subindex;
  
    //return "Overlapping genes between "+ "Cluster " +( d.source.index+1) +" and " +"Cluster " +  (d.source.subindex+1) + ", num of genes:"+dic[index].split("|").length;

});*/


var range = document.getElementById('slider');

noUiSlider.create(range, {
  start: [ maxNumOfOverlappingGens+1 ], // Handle start position
  step: 1, // Slider moves in increments of '10'
  margin: 0, // Handles must be more than '20' apart
  direction: 'rtl', // Put '0' at the bottom of the slider
  orientation: 'horizontal', // Orient the slider vertically
  //behaviour: 'tap-drag', // Move handle on tap, bar is draggable
  range: { // Slider can select '0' to '100'
    'min': minNumOfOverlappingGens,
    'max': maxNumOfOverlappingGens+1
  }
});

var level_g = 0;

var clusterNodesIndex = [];

var inputFormat = document.getElementById('input_slider');

range.noUiSlider.on('update',function(values,handle){
    inputFormat.value = values[handle];
});

range.noUiSlider.on('change', function( values, handle ) {
    var level = getLevelFromNumOfOverlappingGenes(inputFormat.value);
    getClusterNodesIndexBeingSelected(level);

    //temporary fix
    if(labelOff==0)
      updateLabel();
});

/*range.noUiSlider.on('change', function(values) {
    var level = getLevelFromNumOfOverlappingGenes(values);
    getClusterNodesIndexBeingSelected(level);
});*/

inputFormat.addEventListener('change',function(){
    range.noUiSlider.set(this.value);
    var level = getLevelFromNumOfOverlappingGenes(this.value);
    getClusterNodesIndexBeingSelected(level);

    if(labelOff==0)
      updateLabel();
});

var groupText = svg.append("svg:g")
     .selectAll("g")
       .data(chord.groups)
      .enter().append("svg:g");

var text = groupText
    .selectAll("g")
 .data(chord.groups)
    .enter().append("svg:g")
 .attr("transform", function(d) {
   return "rotate(" + ((d.startAngle+d.endAngle)/2 * 180 / Math.PI - 90) + ")"
       + "translate(" + r1 + ",0)";
}).append("svg:text")
   .attr("x", 8)
   .attr("dy", ".45em")
   .attr("text-anchor", function(d) {
     return (d.startAngle+d.endAngle)/2 > Math.PI ? "end" : null;
   })
   .attr("transform", function(d) {
     return (d.startAngle+d.endAngle)/2 > Math.PI ? "rotate(180)translate(-16)" : null;
   })
   .text(function(d) {return (typeof go_inf[d.index].GO_id == "string")? go_inf[d.index].GO_name:getMaxLabel(go_inf[d.index].GO_name)+"+"});



drawLable();


$('.radioButton').change(function(){
  switch(this.value) {
        case "0" :
            drawLable();
            break;
        case "1" :
            drawHierCluster();
            break;
    }            
});

$('#filter').keydown(function(event ){
  if(event.which == 13)
    refreshDetailPanel();
});



$('#filter_button').click(function(){
  refreshDetailPanel();
});

var details_opened = true;

$('#arrow').click(function(){
  toggleDetails();

});

function toggleDetails(){
  $('#arrow').css('transform', function(){ return details_opened ? 'rotate(0deg)' : 'rotate(180deg)'})
  
  $('#details').css('margin-right', function(){ return details_opened ? '-475px' : '0'});
  details_opened = !details_opened;
}
