(function(){

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

  Array.prototype.unique = function() {
      var o = {}, i, l = this.length, r = [];
      for(i=0; i<l;i+=1) o[this[i]] = this[i];
      for(i in o) r.push(o[i]);
      return r;
  };

  var MonaGO = function(){
    var force;
    var nodesSize = size;
    var maxNodeIndex = size;

    var old_array_order = [];
    var clusterHierDataFiltered = [];
    var clusterHierDataStatic = [];

    var memCache = {};
    

    var labelOff=0;

    var clusterHierNodesStatus = [];//for storing nodes generated Hierarchical clustering
    var chordGroupsNodePosition = [];

    var maxNumOfOverlappingGens = clusterHierData[0][3];
    var minNumOfOverlappingGens = clusterHierData[clusterHierData.length-1][3];

    this.gene_info = {};//gene name for each GO term e.g. gene_info[i] = [...]
    this.dic = {};// gene intersection information
    this.go_inf = [];
    var that = this;

    var maxpVal;
    var minpVal;
    var pVal_sort=[];
    var range;
    var step;
    var seperated_points = [];
    var panel=[];

    var fill = d3.scale.ordinal()
         .domain(d3.range(12).reverse())
         .rangeRoundPoints([1,255]);

    var color_scheme = ["#1249C9","#0F399B","#0A2B76","#0F2147","#2C3645","#82733D","#E1C03B","#E6AE29","#F2AB1C","#FFAB00"].reverse();
    var pValLevel = ["p-1","p-2","p-3","p-4","p-5","p-6","p-7","p-8","p-9","p-10"];
    var preElement = "p-1";

    var chordMatrix;
    var chord;
    var w = $(window).width()-475,
         h = $(window).height(),
         r0 = Math.min(w, h) * .25,
         r1 = r0 * 1.1;



    var svg ;
    var go_chart;
    var circleSvg;
    var goLabelSize = {};

    var width = 400;
    var height = 600;
    var force = d3.layout.force()
      .charge(-3000)
      .linkDistance(50)
      .theta(0.2)
      .gravity(1)
      .size([width, height]);

    var zoom = d3.behavior.zoom().translate([w/2, h/2]);
    var zoomLevel = 1;

    
    var level_g = 0;
    var clusterNodesIndex = [];
    var details_opened = true;

    var clusterArc;
    var clusterLine1;
    var clusterLine2;
    var clusterHierNodeView;
    var clusterHierNodeTextView;
    var groupElements;
    var groupLayout;
    var groupText;
    var chordElements;
    var chordLayout;
    var textBackground;


    function createGeneInfo(){
      gene_info = {};
      that.go_inf.forEach(function(d,i){
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

      return gene_info;
    }

    function createInteractionGenes(){
       dic = {};
       var size = that.go_inf.length;
       for(var i = 0;i < size;i++){
          for(var j = 0;j < size;j++){
            if(i == j){
              dic[i+"-"+j] = ""
            }else{
              var genes_i = that.go_inf[i]["genes"];
              var genes_j = that.go_inf[j]["genes"];
              var intersectionGenes = genes_i.intersection(genes_j);
              

              dic[i+"-"+j] = intersectionGenes.getGeneString("|");
            }
          }
       }
       return dic
    }

    function GOFilter(){
      var go_filter = [];
      that.go_inf.map(function(d,i){
          if(Number(d.pVal) < 0.05){
            go_filter.push(d);//remove
          }
      });
      that.go_inf = go_filter;
      nodesSize = that.go_inf.length;
    };

    function getMinPval(){
      var min = Number(that.go_inf[0].pVal);
      for(i=0;i<nodesSize;i++){
        if(Number(that.go_inf[i].pVal) < min){
          min = Number(that.go_inf[i].pVal);

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
        });

       
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

      data = createD3Structure(struct,goid);

      go_chart.selectAll("line.link").remove();
      go_chart.selectAll("circle.node").remove();
      go_chart.selectAll("text.node").remove();

      drawGraph(data);
    }

    function getMaxPval(){
      var max = Number(that.go_inf[0].pVal);
      for(i=0;i<nodesSize;i++){
        if(Number(that.go_inf[i].pVal) > max){
          max = Number(that.go_inf[i].pVal);
        }
      }
      return max;
    }

    function sortGO_inf(){

      pVal_sort = [];
      for(i=0;i<nodesSize;i++)
          pVal_sort[i] = Number(that.go_inf[i].pVal);

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

      return pVal_sort;

    }

    function determineLabelSize(){
      var GOMinSize = 15;
      var GOMaxSize = 30;

      that.go_inf.map(function(d){
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

        for(i=0;i<seperated_points.length-1;i++){
          if((seperated_points[i] <= d)&&(d < seperated_points[i+1])){
              return color_scheme[i];
          }
        }
        if(seperated_points[seperated_points.length-1] == d)
          return color_scheme[9];
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
        panel.push('<i id="'+pValLevel[i]+'" style="background:' + getColor(from,fill) + '" ></i> ' +
                    scienceFormat(from) +" - "+ scienceFormat(to)
        );
      }


      panelLabel = panel.join('<br>');
      $('#pval-label').append(panelLabel);


      function scienceFormat(value){
        return value.toFixed(4).toString();
      }
    }

    function transformIndex(index){
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
          if(zoom.scale()<0.5){
            //reset zoom scale and do nothing
             zoom.scale(0.5)
          }else{
            (transition ? circleSvg.transition() : circleSvg)
              .attr("transform", "translate(" + zoom.translate() + ") scale(" + zoom.scale() + ")");

            if(textBackground){
                if((zoom.scale()<0.8)&&(0.8<zoomLevel<0.9)){
                  textBackground.attr("visibility","hidden");
            }

                if(zoom.scale()>=0.7){
                      textBackground.attr("visibility","visible");
                }
            }


            zoomLevel = zoom.scale();
          }
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

    function getClusterGenes(nodePositions){
      var clustergenes = [];

      nodePositions.map(function(d){ 
            that.go_inf[d].genes.map(function(d,i){
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
      var removeGO = that.go_inf.splice(nodePositions[0],nodePositions.length);
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

      that.go_inf.splice(nodePositions[0],0,newGONode);
      //console.log(go_inf);

      return removeGO;
    }

    function expandNodeSet(nodeBeingClicked,removeGO){
      var position;

      that.go_inf.map(function(d,i){
          if(d["nodeBeingClicked"] == nodeBeingClicked)
          {
            position = i;
          }
      });
      that.go_inf.splice(position,1);

      removeGO.map(function(d,i){
          that.go_inf.splice(position+i,0,d);
      });
      
    }

    function updateMatrix(){
      var newMatrix = [];
      var newArray = [];

      that.go_inf.map(function(d1,i){
          newArray = [];
          that.go_inf.map(function(d2,j){
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
         .style("fill", function(d) { return getColor(Number(that.go_inf[d.index].pVal),fill); })
         .style("stroke", function(d) { return getColor(Number(that.go_inf[d.index].pVal),fill); })
         .attr("d", d3.svg.arc().innerRadius(r0).outerRadius(r1)) //important
         .attr("id","chordGroups")
         .on("mouseover", mouseover_group);

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

      //transform the index to clusterlevel

      var nodeBeingClicked = d["index"];

      if(d["status"]=="Expanded"){

        var nodes = getHierNodes(nodeBeingClicked);

        var nodePositions = getLeafNodesPosition(nodes["leafNodes"]).unique();

        var clusterNodesLevel = nodes["clusterNodesLevel"];

        var genes = getClusterGenes(nodePositions);

        var removeGO = collapseNodeSet(nodePositions,genes,nodeBeingClicked);

        createInteractionGenes();
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

        createInteractionGenes();

        //createGeneInfo();

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

    function createGeneListHtml(index){
      var geneListInHtml = "";

      geneListInHtml += "<div class='gene_content'>";
      that.go_inf[index].genes.forEach(function(d,i){
           var tmp = "<li>" + (i+1) +"."+"<n class='gene_name'>"+ d + "</n> "+ "</li>";
           geneListInHtml+=tmp;
      });
      geneListInHtml += "</div>";

      return geneListInHtml;
    }

    function getNumOfGOTerms(goid){

        if(typeof goid === "string" ){
          return 1;
        }

        return goid.length;
    }

    function createGOList(go_names){
      var goList = "<div class='go_List'>";
      go_names.forEach(function(d,i){
           var tmp = "<li>" + (i+1) +"."+"<n class='Go_name'>"+ d + "</n> "+ "</li>";
           goList+=tmp;
      });
      goList += "</div>";
      return goList;
    }

    function createDetailPanelTempl(i){
        var detailPanelTempl = "";

        var numOfGOTerms = getNumOfGOTerms(that.go_inf[i].GO_id);
        console.log(numOfGOTerms);
        var goInfTempl = "<p> <a class='prop-field'> GO_id: </a>" + that.go_inf[i].GO_id + "</p>";

        if(numOfGOTerms == 1){
           goInfTempl += "<p> <a class='prop-field'>GO_Name: </a>"+ that.go_inf[i].GO_name + "</p>";
        }else{
           goInfTempl += "<a class='prop-field go_dropmenu'> GO_name: <b id='caret_GO' class='caret rotate180'></b></a>" + 
                          createGOList(that.go_inf[i].GO_name) + "<p>";
        }

         goInfTempl += " <p> <a class='prop-field'>Num of genes: </a>"+ that.go_inf[i].count + "</p>" +
                       "<p> <a class='prop-field'>P-value: </a>" + that.go_inf[i].pVal + "</p>";

        var geneListInHtml = createGeneListHtml(i);
        var genesListTempl = "<a class='prop-field gene_dropmenu'>Genes:<b id='caret_gene' class='caret rotate180'></b></a>"+geneListInHtml+"</p>";


        var chartTempl = (numOfGOTerms == 1)?"<p><a class='prop-field'>GO Hierarchy:</a></p> <div id='go_chart'></div> ":"";

        detailPanelTempl += goInfTempl + genesListTempl + chartTempl;

        return detailPanelTempl;
    }

    function setUpDetailPanelListener(){
        var genes_shown = true;
        var go_shown = true;
        $('.gene_dropmenu').click(function(d){
              
          $('#caret_GO').css('transform', function(){ return genes_shown ? 'rotate(0deg)' : 'rotate(180deg)'})
          genes_shown = !genes_shown;

          $geneList = $(this).next();
          $geneList.slideToggle(500);
        });

        $('.gene_name').click(function(){
            $('#filter').val($(this).html());
            refreshDetailPanel();
        });

        $('.go_dropmenu').click(function(d){
              
          $('#caret_gene').css('transform', function(){ return go_shown ? 'rotate(0deg)' : 'rotate(180deg)'})
          go_shown = !go_shown;

          $goList = $(this).next();
          $goList.slideToggle(500);
        });
    }

    function mouseover_group(d, i) {
        $('#content').empty();
        chordLayout.classed("fade", function(p) {
          return p.source.index != i
              && p.target.index != i;
        });
        
        var detailPanelTempl = createDetailPanelTempl(i);
        $('#content').append(detailPanelTempl);

        setUpDetailPanelListener();

        if(typeof that.go_inf[d.index].GO_id == "string"){
          $('#content').append("<div id=\"go_chart\"></div>");
          createGoHier(that.go_inf[d.index].GO_id);
        }

        changePvalPanel(d.index);
    }

    function changePvalPanel(index){
        var levelElement = getPvalLevel(that.go_inf[index].pVal);
        setlevelElement(levelElement);
    }

    function getPvalLevel(pVal){
        for(i=0;i<seperated_points.length-1;i++){
          if((seperated_points[i] <= pVal)&&(pVal < seperated_points[i+1])){
              return pValLevel[i];
          }
        }
        if(seperated_points[seperated_points.length-1] == pVal)
          return pValLevel[9];
    }

    function setlevelElement(element){

      $('#'+preElement).css("border","0");
      $('#'+element).css("border","solid 1px red");
      preElement = element;

    }

    function mouseover_chord(d, i) {
        $('#content').empty();
        var index = d.source.index+"-"+d.source.subindex;

        var str = "";

        dic[index].split("|").forEach(function(d,i){
           var tmp = "<p>" + (i+1) +"."+ d + "</p>";
           str+=tmp;
        });

        $('#content').append("<p>Overlapping genes between " + that.go_inf[d.source.index].GO_id +" and " +  that.go_inf[d.source.subindex].GO_id+":\n</p>"+str);
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
      var genes = "<div class='gene_content'>";

      that.go_inf[num].genes.forEach(function(d,i){

           var tmp = "<li>" + (i+1) +".<n class='gene_name'>"+ d + "</n></li>";
           genes+=tmp;
        });

      genes += "</div>";

      return genes;
    }

    function renderGOTerm(go_num_array){
      var GOterms="";
      for(i=0;i<go_num_array.length;i++){
        num = go_num_array[i];
        Go = that.go_inf[num].GO_id + " " + that.go_inf[num].GO_name;

        var line = "<button class=\"dropbtn\""+ "id=GO_button_"+ num +">" + Go + "</button>";
        line += "<div class=\"Go_content\">" + "<p><a class='prop-field'> P-value:</a> " + that.go_inf[num].pVal + 
        " <p><a class='prop-field'> Num of genes: </a>"+ that.go_inf[num].count +
        "<p><a class='prop-field'>Genes: </a>" +  getGenesFromACluster(num) + "</div>";
        GOterms += line; 

      }
      return GOterms;
    }

    function getGoIndex(go_id){
      for(i in that.go_inf){
          if(that.go_inf[i].GO_id == go_id)
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
        pa = Number(that.go_inf[a].pVal);
        pb = Number(that.go_inf[b].pVal);
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

        for(var i=0;i<that.go_inf.length;i++){
          if(that.go_inf[i]["GO_name"] == searchTerm.toLowerCase()){
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

          for(var j=0;j<that.go_inf.length;j++){
            if(that.go_inf[j]["GO_name"] == searchTerm.toLowerCase()){
              i = j;
            }
          }
        }
        

        if(i != undefined){

          genes = getGenesFromACluster(i);
        
          var detailPanelTempl = createDetailPanelTempl(i);
          $('#content').append(detailPanelTempl);

          setUpDetailPanelListener();

          createGoHier(that.go_inf[i].GO_id);

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
        var indexs = that.gene_info[searchTerm];

        
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
            $('.gene_name').click(function(){
                $('#filter').val($(this).html());
                refreshDetailPanel();
            });

            createOnClick(num_array);

            var arrayList=[];
            for(i in index_array){
              var num = parseInt(index_array[i]);
              arrayList.push(num);
            }

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
    }

    function updateLabel(){
      groupText.remove();

      groupText = circleSvg.append("svg:g")
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
            .text(function(d) {return (typeof that.go_inf[d.index].GO_id == "string")? that.go_inf[d.index].GO_name:getMaxLabel(that.go_inf[d.index].GO_name)+"+"});
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

    function drawHierCluster(){
      labelOff = 1;
      clusterArc.style("display","");
      clusterLine1.style("display","");
      clusterLine2.style("display","");
      clusterHierNodeView.style("display","");
      clusterHierNodeTextView.style("display","");


      groupText.style("display","none");

      drawHierarchicalClustering(clusterHierData);
      if(textBackground){
        if(zoom.scale()<0.7){
            textBackground.attr("visibility","hidden");
        }else{
            textBackground.attr("visibility","visible");
        }
      }
    }

    function toggleDetails(){
      $('#arrow').css('transform', function(){ return details_opened ? 'rotate(0deg)' : 'rotate(180deg)'})
      
      $('#details').css('margin-right', function(){ return details_opened ? '-475px' : '0'});
      details_opened = !details_opened;

      redrawMain_vis(details_opened);
      movePvalPanel(details_opened);
    }

    function redrawMain_vis(details_opened){

      if(!details_opened){
        d3.select(".main_vis").attr("width",w+475);
        circleSvg.transition().attr("transform", "translate(" + (w+475) / 2 + "," + h / 2 + ") scale(" + zoom.scale() + ")");
      }else{
        d3.select(".main_vis").attr("width",w);
        circleSvg.transition().attr("transform", "translate(" + w / 2 + "," + h / 2 + ") scale(" + zoom.scale() + ")");
      }
    }

    function movePvalPanel(details_opened){
      if(!details_opened){
        d3.select(".pval-label").transition().style("margin-left",-190);
      }else{
        d3.select(".pval-label").transition().style("margin-left",-660);
      }
    }

    function setUpControlPanel(){

      var element = '<label style="font-size:16px">Cluster GO term according to the minimum number of common gene(s)</label> \
          <div id="slider" class="sliderBar"></div>\
          <input type="text" id="input_slider"/>\
            <table class="RadioBox">\
            <tbody><tr><td>\
                  <input id="labelRadioBox" class="radioButton" type="radio" name="radioBox" value="0" checked>\
                </td><td><label for="labelRadioBox">Show GO term name</label>\
                </td></tr><tr><td></td></tr><tr><td>\
                  <input id="hierClusterRadioBox"  class="radioButton" type="radio" name="radioBox" value="1" >\
                </td><td>\
                   <label style="width:300px;" for="hierClusterRadioBox">Show hierarchical tree and click on the node to manually cluster the GO term</label>\
                </td></tr></tbody></table>';

      $("#control-panel").append(element);
      var ranger = document.getElementById('slider');
      var inputFormat = document.getElementById('input_slider');
      setUpRangeSlider(ranger,inputFormat,minNumOfOverlappingGens,maxNumOfOverlappingGens);
    }

    function setUpRangeSlider(ranger,inputFormat,minNumOfOverlappingGens,maxNumOfOverlappingGens){
      noUiSlider.create(ranger, {
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

      ranger.noUiSlider.on('update',function(values,handle){
          inputFormat.value = values[handle];
      });

      ranger.noUiSlider.on('change', function( values, handle ) {
          var level = getLevelFromNumOfOverlappingGenes(inputFormat.value);
          getClusterNodesIndexBeingSelected(level);

          //temporary fix
          if(labelOff==0)
            updateLabel();
      });

      inputFormat.addEventListener('change',function(){
          ranger.noUiSlider.set(this.value);
          var level = getLevelFromNumOfOverlappingGenes(this.value);
          getClusterNodesIndexBeingSelected(level);

          if(labelOff==0)
            updateLabel();
      });
    }

    function setUpListener(){
      $('#filter').keydown(function(event){
        if(event.which == 13){
          refreshDetailPanel();
          $('#searchBox').remove();
        }

      });

      $('#filter_button').click(function(){
        refreshDetailPanel();
        $('#searchBox').remove();
      });

      $('#arrow').click(function(){
        toggleDetails();
      });

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
    }

    function setUpView(){
      svg = d3.select("#chart")
       .append("svg:svg")
         .attr("class","main_vis")
         .attr("width", w)
         .attr("height", h)
         .call(zoom.on("zoom", redraw))

      circleSvg = svg.append("svg:g")
         .attr("id", "circle")
         .attr("transform", "translate(" + w / 2 + "," + h / 2 + ")")

      circleSvg.append("circle").attr("r",r0);
      
      clusterArc = circleSvg.append("svg:g")
         .selectAll("g")
           .data(clusterHierData)
         .enter().append("path");

      clusterLine1 = circleSvg.append("svg:g")
       .selectAll("g")
         .data(clusterHierData)
       .enter().append("path");

      clusterLine2 = circleSvg.append("svg:g")
       .selectAll("g")
         .data(clusterHierData)
       .enter().append("path");

      clusterHierNodeView = circleSvg.append("svg:g")
       .selectAll("g")
         .data(clusterHierNodesStatus)
       .enter().append("svg:g")
        .attr("id","clusterNode")

      clusterHierNodeTextView = circleSvg.append("svg:g")
       .selectAll("g")
        .data(clusterHierNodesStatus)
       .enter().append("svg:g");

      groupElements = circleSvg.append("svg:g")
       .selectAll("path")
         .data(chord.groups)
         .enter().append("svg:path");

      groupLayout = groupElements
         .style("fill", function(d) { return getColor(Number(that.go_inf[d.index].pVal),fill); })
         .style("stroke", function(d) { return getColor(Number(that.go_inf[d.index].pVal),fill); })
         .attr("d", d3.svg.arc().innerRadius(r0).outerRadius(r1)) //important
         .attr("id","chordGroups")
         .on("mouseover", mouseover_group);

      groupText = circleSvg.append("svg:g")
           .selectAll("g")
             .data(chord.groups)
            .enter().append("svg:g");

      that.go_inf.map(function(d){
        if(goLabelSize[d.GO_name] == undefined){
          goLabelSize[d.GO_name] = d.size;
        }
      })

      chordElements = circleSvg.append("svg:g")
            .selectAll("path")
            .data(chord.chords)
            .enter().append("svg:path");

      chordLayout = chordElements
            .attr("class", "chord")
            //.style("fill", function(d) { return fill(d.source.index%10); })
            //.style("opacity",".5")
            .attr("d", d3.svg.chord().radius(r0))
            .attr("id","chordChords")
            .on("mouseover", mouseover_chord);

      setUpControlPanel();
    }

    this.init = function(size,go_inf,clusterHierData){

      that.go_inf = go_inf;

      clusterHierData.map(function(d){
        clusterHierDataStatic.push([d[2],d[3]]);
      });

      that.go_inf.forEach(function(d,i){
        that.go_inf[i].genes = d.genes.split("|");
      });


      maxpVal = getMaxPval();
      minpVal = getMinPval();

      pVal_sort = sortGO_inf();

      step = (maxpVal-minpVal)/10;

      chordMatrix = d3.layout.chord()
       .padding(.03)
       .sortSubgroups(d3.descending);
       
      chord = chordMatrix.matrix(matrix);

      //store the nodes for hierarchical clustering visualzaition
      chord.groups().forEach(function(d,i){
        nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
        chordGroupsNodePosition.push(nodeObj);
      });
       
      for(i=0;i<11;i++){
        seperated_points.push(minpVal+step*i);
      }



      that.gene_info = createGeneInfo();
      that.dic = createInteractionGenes();

      createPvalLabelPanel(fill);
      createClusterNodesStatus(clusterHierData,[],"",[]);
      setUpView();
      determineLabelSize();
      drawLable();
      setUpListener();
      
    }

  }

  new MonaGO().init(size,go_inf,clusterHierData);


})();











