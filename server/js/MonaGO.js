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

    var monago;

    var force;
    var nodesSize = size;
    var maxNodeIndex = size;

    var old_array_order = [];
    var clusterHierDataFiltered = [];

    
    

    var labelOff=0;

    var clusterHierNodesStatus = [];//for storing nodes generated Hierarchical clustering
    var chordGroupsNodePosition = [];



    var popUpList = [];//store the index of go that being popped up.

    this.gene_info = {};//gene name for each GO term e.g. gene_info[i] = [...]
    this.dic = {};// gene intersection information
    this.go_inf = [];
    this.goNodes = [];
    this.matrix = [];
    this.clusterHierData = [];
    this.clusterHierDataStatic = [];
    this.memCache = {};
    this.level_g = 0;

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
    var detailPanelWidth = $(window).width()*0.25;

    var controlPanelWidth = 520;

    var controlPanelHeight = 120;

    var w = $(window).width()-detailPanelWidth,
         h = $(window).height(),
         r0 = Math.min(w, h) * .20,
         r1 = r0 * 1.1;

    var goNameArr = [];


    var main_div = 0;
    var svg ;
    var go_chart;
    var circleSvg;
    var goLabelSize = {};

    var width = detailPanelWidth-70;
    var height = 600;
    var force = d3.layout.force()
      .charge(-3000)
      .linkDistance(50)
      .theta(0.2)
      .gravity(1)
      .size([width, height]);

    var zoom = d3.behavior.zoom().translate([w/2, h/2]);
    var zoomLevel = 1;

  
    var clusterNodesIndex = [];
    var details_opened = true;
    var control_opened = true;
    var reloadData = false;

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
    var app = angular.module('MonaGO');

    app.directive('clearBtn', ['$parse', function ($parse) {
        return {
            link: function (scope, elm, attr, ngModelCtrl) {

                elm.wrap("<div style='position: relative'></div>");
                var btn = '<i class="searchclear ng-hide fa fa-close"></i>';
                var angularBtn = angular.element(btn);

                angularBtn.css('top', 15);
                elm.after(angularBtn);

                //clear the input
                angularBtn.on("click", function () {
                    popUpList = [];
                    scope.searchText = "";
                    $('#filter').val("");
                    $('#filter').focus();
                    refreshDetailPanel();
                    if($('#searchBox')) $('#searchBox').remove();
                    angularBtn.addClass("ng-hide");
                    scope.$apply();
                });

                // show  clear btn  on focus
                elm.bind('focus keyup change paste propertychange', function (blurEvent) {
                    if (elm.val() && elm.val().length > 0) {
                        angularBtn.removeClass("ng-hide");

                    } else {
                        angularBtn.addClass("ng-hide");
                    }
                });

            }
        };
    }]);


    //associate genes with go,cluster(using index in go_inf)
    function createGeneInfo(){
      that.gene_info = {};
      that.go_inf.forEach(function(d,i){
        var index = i;
        d["genes"].forEach(function(d,i){
          if(that.gene_info[d]!=undefined){
              var str =   ";"  +  index.toString();
              that.gene_info[d] = that.gene_info[d] + str;
          }
          else{
            that.gene_info[d] = index.toString();
          }
        });
      });

    }

    function createInteractionGenes(){
       that.dic = {};
       var size = that.go_inf.length;
       for(var i = 0;i < size;i++){
          for(var j = 0;j < size;j++){
            if(i == j){
              that.dic[i+"-"+j] = ""
            }else{
              var genes_i = that.go_inf[i]["genes"];
              var genes_j = that.go_inf[j]["genes"];
              var intersectionGenes = genes_i.intersection(genes_j);
              
              that.dic[i+"-"+j] = intersectionGenes.getGeneString("|");
            }
          }
       }
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
      return that.goNodes[goid];
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
    function drawGraph(struct,svg) {
      force
        .nodes(struct.nodes)
        .links(struct.links);

      var strokeColor = d3.scale.category10();
      var radiusScale = d3.scale.linear().domain([0,d3.max(struct.nodes,function(d) {return d.r;})]).range([5,10]);

      var link = svg.selectAll("line.link")
        .data(struct.links)
        .enter().append("line")
        .attr("class", "link")
        .style('stroke', function(d, i ) { return strokeColor(d.type);})
        .style("stroke-width", 2);
     
      var node = svg.selectAll("circle.node")
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

       
      var nodeText = svg.selectAll("text.node")
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
    function update(goid,svg,width,height) {

      force = {};
      
      force = d3.layout.force()
      .charge(-3000)
      .linkDistance(50)
      .theta(0.2)
      .gravity(1)
      .size([width, height]);

      struct = parentGO(goid);

      data = createD3Structure(struct,goid);

      svg.selectAll("line.link").remove();
      svg.selectAll("circle.node").remove();
      svg.selectAll("text.node").remove();

      drawGraph(data,svg);
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
      var panelTitle = [];
      var panel = [];
      panelTitle.push("p-value");
      panelTitle.push("<hr></hr>");
      //create panel
      for(i=0;i<seperated_points.length-1;i++)
      {
        from = seperated_points[i];
        to = seperated_points[i+1];
        panel.push('<i id="'+pValLevel[i]+'" style="background:' + getColor(from,fill) + '" ></i> ' +
                    scienceFormat(from) +" - "+ scienceFormat(to)
        );
      }


      panelLabel = "<div>" + panelTitle + panel.join('<br>') + "</div>";
      $('#pval-label').empty();
      $('#pval-label').append(panelLabel);


      function scienceFormat(value){
        return value.toFixed(4).toString();
      }
    }

    function transformIndex(index){

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

    function drawClusterNodes(){
      var clusterNodes = [];
      clusterHierNodesStatus.map(function(d){
        if(d["numOfOverlappingGenes"]>0)
          clusterNodes.push(d);
      })

      clusterHierNodeView = clusterHierNodeView.data(clusterHierNodesStatus,function(d){return d.index;});

      clusterHierNodeView.enter().append("g");

      clusterHierNodeView.exit().remove();

      clusterHierNodeView.transition().duration(500).attrTween("transform", tween).attr("class","clusterNodeView");

      clusterHierNodeView.selectAll("circle").remove();

      circle = clusterHierNodeView.append("circle");

      circle.attr("r",function(d){
        if((d["numOfOverlappingGenes"] > 0)||(d["numOfOverlappingGenes"]==undefined))
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
          if((d["numOfOverlappingGenes"] > 0)||(d["numOfOverlappingGenes"]==undefined))
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
          if((d!=undefined) && (d[3]!=0))
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
          if((d!=undefined) && (d[3]!=0))
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
        drawClusterNodes();
        
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

   /*get all the leaf nodes along with thier cluster levels
    *@parama index of clicked cluster node
    *@return index of leaf nodes and their cluster levels
    */ 
    function getHierNodes(index){

      var leafNodes = [];
      var clusterNodesLevel = [];

      recursiveGetNodes(index-maxNodeIndex);

      function recursiveGetNodes(pos){

          var leftNode = that.clusterHierData[pos][0];
          var rightNode = that.clusterHierData[pos][1];



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


    /**
    *collapse/combine all go inf assoicated with the clustered node and return the removed go inf
    *@parma leaf node position in chord(get all the assoicated go inf),total genes in the cluster, index of the node being clicked
    *@return original go inf(before clustered)
    */
    function collapseNodeSet(nodePositions,nodeBeingClicked){
      //get all the genes associated with the cluster
      var genes = getClusterGenes(nodePositions);

      var clusterGOs = that.go_inf.splice(nodePositions[0],nodePositions.length);
      var GO_id = [];
      var GO_name = [];
      clusterGOs.map(function(d){
          GO_id.push(d["GO_id"]);
          GO_name.push(d["GO_name"]);
      });

      var pVal = calculateNewPvalue(clusterGOs);

      var newGONode = {"GO_id":GO_id,"GO_name":GO_name,"count":genes.length,"pVal":pVal,"genes":genes,"nodeBeingClicked":nodeBeingClicked};

      that.go_inf.splice(nodePositions[0],0,newGONode);

      return clusterGOs;
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

    function drawLayout(newMatrix){


        groupElements = groupElements.data(chord.groups());

        groupElements.enter().append("svg:path");

        groupElements.exit().remove();

        groupLayout = groupElements
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

    /**
    *create cluster nodes(position,status)
    *@param clusterHierData: hierarchical data
    *@param nodeBeingClicked: index of the node being clicked
    *@param status: expanded/collapsed
    *@param collpaseNodes: an array of collapsed nodes 
    *@param clusterNodesRadius: an array of collapsed nodes' radius
    *@param collapsedNodeRadius: the radius of the collapsed node
    */
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

    /**
    *remove associated cluster in cluster hierarchical data and create new cluster hierarchical data
    *@param nodePositions : the node position assoicated with the collapsed cluster
    *@param clusterNodesLevel : 
    *@return new created cluster hierarchical data
    */
    function createNewClusterHierData(nodePositions,clusterNodesLevel){
        var newNode = old_array_order[nodePositions[nodePositions.length-1]];
        var newClusterHierData = [];

        that.clusterHierData.map(function(d,i){
            if(clusterNodesLevel.indexOf(i) == -1){
                newClusterHierData.push(d);
            }else{//create new node
                if(clusterNodesLevel.indexOf(i) == clusterNodesLevel.length-1)
                {
                    newClusterHierData.push([newNode,newNode,d[2]]);
                }else{
                    newClusterHierData.push(undefined);
                }
            }
        });    

        return newClusterHierData;
    }

    function addNodesInClusterData(removedClusterHierData){

      removedClusterHierData.map(function(d,i){

          that.clusterHierData[d[2]-maxNodeIndex] = d;
      });


      return that.clusterHierData;
    }

    function getRemoveClusterData(clusterNodesLevel){
      var removeHierData = [];

      clusterNodesLevel.map(function(d,i){
          removeHierData.push(that.clusterHierData[d]);
      });

      return removeHierData;
    }

    /**
    *collpase the hierarcial clusters and return removed hierarcial data
    *@param nodePositions: an array of nodes' positions
    *@param clusterNodesLevel: the level of the nodes being clustered in hierarchical data
    *@param nodeBeingClicked: the index of the node being clicked
    *@param collapsedNodes: an array of collapsed nodes
    *@param clusterNodesRadius: an array of collapsed nodes' radius  
    *@return removed hierarcial data
    */
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


        //create new cluster hierarchical cluster data
        that.clusterHierData = createNewClusterHierData(nodePositions,clusterNodesLevel);

        //create cluster node 
        createClusterNodesStatus(that.clusterHierData,nodeBeingClicked,"collapse",collpasedNodes,clusterNodesRadius,collapsedNodeRadius);

        return [removeHierData,that.clusterHierData];
    }

    function expandHierClustering(removedClusterHierData,clusterNodesLevel,nodeBeingClicked,collpasedNodes,clusterNodesRadius){
        clusterHierNodesStatus = [];//empty clusterHierNodesStatus
        chordGroupsNodePosition = [];

        chord.groups().forEach(function(d,i){//recalculate angel,radius
          nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
          chordGroupsNodePosition.push(nodeObj);
        });

        that.clusterHierData = addNodesInClusterData(removedClusterHierData,clusterNodesLevel);

        createClusterNodesStatus(that.clusterHierData,nodeBeingClicked,"expand",collpasedNodes,clusterNodesRadius);

        return that.clusterHierData;
    }

    function collapseMonaGO(nodeBeingClickedIndex){
      //get leaf nodes(associated with the cluster) info(node index & cluster level)
        var nodes = getHierNodes(nodeBeingClickedIndex);

        //get the leaf nodes position in the chord layout
        var nodePositions = getLeafNodesPosition(nodes["leafNodes"]).unique();
        var clusterNodesLevel = nodes["clusterNodesLevel"];

        //key function: collpase/combine go and return the go inf being removed
        var removeGOs = collapseNodeSet(nodePositions,nodeBeingClickedIndex);

        var removeNodeInArray = removeNodesInArrayOrder(nodePositions);

        matrix = updateMatrix();
        chord = chordMatrix.matrix(matrix);

        var collapseResults = collapseHierClustering(nodePositions,clusterNodesLevel,nodeBeingClickedIndex,Object.keys(that.memCache),that.memCache["clusterNodesRadius"]);
        removedHierData = collapseResults[0];
        clusterHierData = collapseResults[1];
        

        if(that.memCache["clusterNodesRadius"]!=undefined)
            var clusterNodesRadius = that.memCache["clusterNodesRadius"];
        else
            clusterNodesRadius = [];

        clusterHierNodesStatus.map(function(d){
            if(d["index"]==nodeBeingClickedIndex){
              //console.log(d);
              clusterNodesRadius.push({"nodeId":d["index"],"radius":d["radius"]});
            }
        });

        that.memCache["clusterNodesRadius"] = clusterNodesRadius;

        //console.log(clusterNodesRadius);
        /*console.log(removeHierData);*/

        var nodeBeingMemorized = {"go_inf":removeGOs,"array_order":removeNodeInArray,"clusterHierData":removedHierData,"clusterNodesLevel":clusterNodesLevel};

        that.memCache[nodeBeingClickedIndex] = nodeBeingMemorized;

        return clusterHierData;

    }

    function expandMonaGO(nodeBeingClickedIndex){

        nodeBeingMemorized = that.memCache[nodeBeingClickedIndex];

        expandNodeSet(nodeBeingClickedIndex,nodeBeingMemorized["go_inf"]);

        addNodesToArrayOrder(nodeBeingMemorized["array_order"]);

        var collapsedNodes = getMemKey(that.memCache);

        var clusterNodesLevel = nodeBeingMemorized["clusterNodesLevel"];

        var clusterNodesRadius = that.memCache["clusterNodesRadius"];

        var newclusterNodesRadius = [];

        clusterNodesRadius.map(function(d,i){
          if(d["nodeId"]!=nodeBeingClickedIndex)
            newclusterNodesRadius.push(d);
        })

        clusterNodesRadius = newclusterNodesRadius;

        that.memCache["clusterNodesRadius"] = clusterNodesRadius;

        matrix = updateMatrix();
        chord = chordMatrix.matrix(matrix);


        that.clusterHierData = expandHierClustering(nodeBeingMemorized["clusterHierData"],clusterNodesLevel,nodeBeingClickedIndex,collapsedNodes,clusterNodesRadius);
        
        delete that.memCache[nodeBeingClickedIndex];

        return that.clusterHierData;
    }

    function onClusterNodeClick(d,i){

      //transform the index to clusterlevel

      var nodeBeingClickedIndex = d["index"];

      if(d["status"]=="Expanded"){

        clusterHierData = collapseMonaGO(nodeBeingClickedIndex);

      }else{
        /*console.log("node going to expand");*/

        clusterHierData = expandMonaGO(nodeBeingClickedIndex);

      }

      updateMoanaGOLayout(clusterHierData);
      updateLabel();
    }

    function updateMoanaGOLayout(clusterHierData){

      createGOAndUpdateChordLayout();

      //if(labelOff==1)
          drawHierarchicalClustering(clusterHierData);
    }

    function createGOAndUpdateChordLayout(){
        createInteractionGenes();
        createGeneInfo();

        drawLayout(that.matrix);
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
      that.clusterHierDataStatic.map(function(d,i){
          if(d[1]>=numOfOverlappingGenes)
            level = i;
      });

      return level;
    }

    function getClusterNodesIndexBeingSelected(level){
      console.log(level);
      console.log(that.level_g);
      var clusterDataLevel = [];
      if(level>=that.level_g){//collapse
        
        for(var i=level;i >= 0; i--){
          if(that.clusterHierData[i]!=undefined&&clusterDataLevel.indexOf(i)==-1){
            if(that.memCache[clusterHierData[i][2]]==undefined){
              var clusterNodesLevel = getHierNodes(that.clusterHierData[i][2])["clusterNodesLevel"];
              //console.log(clusterNodesLevel);
              clusterNodesLevel.map(function(d){return clusterDataLevel.push(d);});
              clusterDataLevel = clusterDataLevel.unique();
              clusterNodesIndex.push(that.clusterHierData[clusterNodesLevel[clusterNodesLevel.length-1]][2]);
            }
          }
          clusterNodesIndex = clusterNodesIndex.unique();


        }

      }else{//expand
        console.log("expand");
        var index = level + nodesSize;
        var pos = [];

        clusterNodesIndex.reverse().map(function(d1,i){
          if(d1 >= index){
            pos.push(i);

            var clusterNode = getClusterNode(d1);
            that.clusterHierData = expandMonaGO(clusterNode["index"]);
          
          }
        });

        clusterNodesIndex.splice(pos[0],pos.length);

        for(var i=level;i >= 0; i--){
          if(that.clusterHierData[i]!=undefined&&clusterDataLevel.indexOf(i)==-1){

            var clusterNodesLevel = getHierNodes(that.clusterHierData[i][2])["clusterNodesLevel"];
            //console.log(clusterNodesLevel);
            clusterNodesLevel.map(function(d){return clusterDataLevel.push(d);});
            clusterDataLevel = clusterDataLevel.unique();
            clusterNodesIndex.push(that.clusterHierData[clusterNodesLevel[clusterNodesLevel.length-1]][2]);

          }
        }
        clusterNodesIndex = clusterNodesIndex.unique();

      }

      if(clusterNodesIndex.length!=0){
        clusterNodesIndex.map(function(index,i){
            clusterHierNodesStatus.map(function(d){
              if(index == d["index"]){
                if(that.memCache[index]==undefined){
                  that.clusterHierData = collapseMonaGO(index);

                }
              }
            });
            
        });
      }



      that.level_g = level;
      updateMoanaGOLayout(that.clusterHierData);
    }

    function getClusterNode(index){
      var node;
      clusterHierNodesStatus.map(function(d){
        if(d["index"] == index)
          node = d;
      });
      return node;
    }

    function createGoHierifNecessary(goid){

      if(typeof goid == "string"){
          $('#content').append("<div id=\"go_chart\"></div>");
      
      var go_chart = d3.select("#go_chart").append("svg")
        .attr("width", width)
        .attr("height", height);
      go_chart.append('rect')
        .style('fill','white')
        .style('stroke','gray')
        .attr('width',width)
        .attr('height',height)
        .attr('x',0)
        .attr('y',0);

      update(goid,go_chart,width,height);

      }

      
    }

    function replaceCommaWithUnderline(term){
        return term.replace(":","_");
    }

    function createGoHierByGOId(goid){

      var goid_target = replaceCommaWithUnderline(goid);

      var go_chart = d3.select('#'+goid_target).append("svg")
        .attr("width", width)
        .attr("height", height);

      go_chart.append('rect')
        .style('fill','white')
        .style('stroke','gray')
        .attr('width',width)
        .attr('height',height)
        .attr('x',0)
        .attr('y',0);

      update(goid,go_chart,width,height);
    }

    function createGeneListHtml(go){
      var geneListInHtml = "";

      geneListInHtml += "<div class='gene_content'>";
      go.genes.forEach(function(d,i){
           var tmp = "<p>" + (i+1) +"."+"<n class='gene_name'>"+ d + "</n> "+ "</p>";
           geneListInHtml+=tmp;
      });
      geneListInHtml += "</div>";

      return geneListInHtml;
    }

    function getGODetailByName(goName){
        var goDetailArr={};

        that.go_inf_ori.forEach(function(d,i){
          if(d.GO_name == goName){
            goDetailArr = d;
          }
        });

        return goDetailArr;
    }

    function createGOdetailTempl(goName){
      var goDetail = "";

      var goDetailArr = getGODetailByName(goName);

      goDetail += "<div class='go_detail_content'>";
      goDetail += "<p> <a class='prop-field'> GO_id: </a>" + goDetailArr.GO_id + "</p>";
      goDetail += " <p> <a class='prop-field'>Num of genes: </a>"+ goDetailArr.count + "</p>" +
                       "<p> <a class='prop-field'>P-value: </a>" + goDetailArr.pVal + "</p>";

      var geneListInHtml = createGeneListHtml(goDetailArr);
      var genesListTempl = "<a class='prop-field gene_dropmenu'>Genes:</a><b id='caret_gene' class='caret rotate180'></b>"+geneListInHtml+"</p>";

      goDetail += genesListTempl;
      goDetail += "<div id='"+ replaceCommaWithUnderline(goDetailArr.GO_id) +"'>" +"</div>";

      goDetail += "</div>";

      return goDetail;
    }

    function getNumOfGOTerms(goid){

        if(typeof goid === "string" ){
          return 1;
        }

        return goid.length;
    }

    function createGOList(go_names){
      var goList = "<div class='go_List'>";
      goNameArr = [];

      go_names.forEach(function(d,i){
           var goDetail = createGOdetailTempl(d);

           goNameArr.push(d);

           var tmp = "<li>" + (i+1) +"."+"<n class='Go_name go_detail_dropmenu'>"+ d + 
           "</n> <b id='caret_GO_details' class='caret'></b>"+ goDetail + "</li>";
           goList += tmp;
      });


      goList += "</div>";
      return goList;
    }

    function recursiveGetArrays(arrsSource){
      var arrDest = [];

      recursiveGet(arrsSource);

      function recursiveGet(arrsSource){
        if (typeof arrsSource == 'string'){
          arrDest.push(arrsSource);
        }else{
          arrsSource.forEach(function(d){
              recursiveGet(d);
          });
        }
      }
      
      return arrDest;
    }

    function createDetailPanelTempl(i){
        var detailPanelTempl = "";

        var numOfGOTerms = getNumOfGOTerms(recursiveGetArrays(that.go_inf[i].GO_id));
        var goInfTempl = "<p> <a class='prop-field'> GO_id: </a>" + recursiveGetArrays(that.go_inf[i].GO_id).join(", ") + "</p>";

        if(numOfGOTerms == 1){
           goInfTempl += "<p> <a class='prop-field'>GO_Name: </a>"+ recursiveGetArrays(that.go_inf[i].GO_name) + "</p>";
           goInfTempl += " <p> <a class='prop-field'>Num of genes: </a>"+ that.go_inf[i].count + "</p>"+
           "<p> <a class='prop-field'>P-value: </a>" + that.go_inf[i].pVal + "</p>";
        }else{
           goInfTempl += "<a class='prop-field go_dropmenu'> GO_name: <b id='caret_GO' class='caret rotate180'></b></a>" + 
                          createGOList(recursiveGetArrays(that.go_inf[i].GO_name)) + "<p>";
           goInfTempl += " <p> <a class='prop-field'>Num of genes: </a>"+ that.go_inf[i].count + "</p>"+
           "<p> <a class='prop-field'>P-value(Average): </a>" + that.go_inf[i].pVal + "</p>";
        }
           

        var geneListInHtml = createGeneListHtml(that.go_inf[i]);
        var genesListTempl = "<a class='prop-field gene_dropmenu'>Genes:</a><b id='caret_gene' class='caret rotate180'></b>"+geneListInHtml+"</p>";


        var chartTempl = (numOfGOTerms == 1)?"<p><a class='prop-field'>GO Hierarchy:</a></p> <div id='go_chart'></div> ":"";

        detailPanelTempl += goInfTempl + genesListTempl + chartTempl;

        return detailPanelTempl;
    }

    function getCaretStatus(element){

      if(d3.select(element).node().attr("status")==undefined){
        return false;
      };

      if(d3.select(element).node().attr("status")=="true"){
        return true;
      }

      return false;
    }

    function setCaretStatus(element,status){
      if (status == true){
        d3.select(element).node().attr("status","false");
      }else{
         d3.select(element).node().attr("status","true");
      }
    }

    function setUpDetailPanelListener(){
        var genes_shown = true;
        var go_shown = true;
        var go_detail_shown = true;

        $('.gene_dropmenu').click(function(d){
          
          var status = getCaretStatus($(this));

          $(this).next().css('transform', function(){ return !status ? 'rotate(0deg)' : 'rotate(180deg)'});
          setCaretStatus($(this),status);

          $geneList = $(this).next().next();
          $geneList.slideToggle(500);
        });

        $('.gene_name').click(function(){
            $('#filter').val($(this).html());
            $('#filter').focus();
            refreshDetailPanel();
        });

        $('.go_id').click(function(){
            $('#filter').val($(this).html());
            refreshDetailPanel();
        });
        

        $('.go_dropmenu').click(function(d){
              
          $('#caret_GO').css('transform', function(){ return go_shown ? 'rotate(0deg)' : 'rotate(180deg)'});
          go_shown = !go_shown;

          $goList = $(this).next();
          $goList.slideToggle(500);
        });

        $('.go_detail_dropmenu').click(function(d){

          var status = getCaretStatus($(this));
          $(this).next().css('transform', function(){ return status ? 'rotate(0deg)' : 'rotate(180deg)'});
          setCaretStatus($(this),status);

          $goList = $(this).next().next();
          $goList.slideToggle(500);
        });
        
    }

    function createGOHierForClusterGO(){
      goNameArr.forEach(function(d,i){

        goDetail = getGODetailByName(d)
        var goNameID = replaceCommaWithUnderline(goDetail.GO_id);

        var hierHeight = height-20;
        var hierWidth = width-50;

        var go_chart = d3.select("#"+goNameID).append("svg")
        .attr("width", hierWidth)
        .attr("height", hierHeight);

        go_chart.append('rect')
          .style('fill','white')
          .style('stroke','gray')
          .attr('width',hierWidth)
          .attr('height',hierHeight)
          .attr('x',0)
          .attr('y',0);

        update(goDetail.GO_id,go_chart,hierWidth,hierHeight);

      })
    }

    function mouseover_group(d, i) {
        $('#content').empty();
        chordLayout.classed("fade", function(p) {
          return p.source.index != i
              && p.target.index != i;
        });

        highlight_index = i;
        groupText.transition().duration(300).style("fill", function(d,i) {
          if(d.index == highlight_index){
            return "black";
          }else{
            return "grey";
          }
        });
        
        var detailPanelTempl = createDetailPanelTempl(i);
        $('#content').append(detailPanelTempl);

        setUpDetailPanelListener();

        createGoHierifNecessary(that.go_inf[i].GO_id);
        createGOHierForClusterGO();

        changePvalPanel(d.index);
    }

    function mouseleave_group(){

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
      $('#'+element).css("border","solid 3px red");
      preElement = element;

    }


    function createChordTempl(index,d){
        
        var geneListInHtml = "<div class='gene_content'>";
        that.dic[index].split("|").forEach(function(d,i){
             var tmp = "<li>" + "<a class='gene_name'>"+ (i+1) +"."+ d + "</a> "+ "</li>";
             geneListInHtml+=tmp;
        });
        geneListInHtml += "</div>";

        templ="<p>Overlapping genes between <a class='go_id'>" + that.go_inf[d.source.index].GO_id + 
        "</a> and <a class='go_id'>" +  that.go_inf[d.source.subindex].GO_id+"</a>:\n</p>"+geneListInHtml;

        return templ;
    }

    function mouseover_chord(d, i) {
        $('#content').empty();

        var index = d.source.index+"-"+d.source.subindex;

        var chordTempl = createChordTempl(index,d);

        $('#content').append(chordTempl);
        setUpDetailPanelListener();
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
          $content.slideToggle(500);
        });

    }

    function createOnHover(num_array){

      $(".dropbtn").mouseover(function (d) {
/*        console.log($(this)[0].id);
        console.log(d.target.id);*/
        var num = getTargetGONum(d.target.id);
        groupLayout.classed("highlight", function(p){
          return p.index==num;
        })

        changePvalPanel(num);
      });

      $(".dropbtn").mouseleave(function(d){
          var num = getTargetGONum(d.target.id);
          groupLayout.classed("highlight",false,function(d){
              return d.index==num;
          });

          resetPvaluePanel(num);
      });

    }

    function resetPvaluePanel(index){
      var levelElement = getPvalLevel(that.go_inf[index].pVal);
      $('#'+levelElement).css("border","0");
    }

    function getTargetGONum(id){
      var id_arr = id.split("_");
      return id_arr[id_arr.length-1];
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
      var num = 0;
      var GO = "";
      for(i=0;i<go_num_array.length;i++){
        var num = go_num_array[i];
        GO = that.go_inf[num].GO_id + " " + that.go_inf[num].GO_name;

        var numOfGOTerms = getNumOfGOTerms(that.go_inf[num].GO_id);

        var line = "<button class=\"dropbtn\""+ "id=GO_button_"+ num +">" + GO + "</button>";
        line += "<div class=\"Go_content\">" + "<p><a class='prop-field'>" + ((numOfGOTerms == 1)?" P-value:</a> ":" P-value(Average):</a>")
         + that.go_inf[num].pVal + " <p><a class='prop-field'> Num of genes: </a>"+ that.go_inf[num].count +
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
        groupText.transition().duration(500).attr("x",8);
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

          createGoHierifNecessary(that.go_inf[i].GO_id);

          chordLayout.classed("fade", function(p) {
                  return p.source.index != i
                      && p.target.index != i;
                });

          groupLayout.transition().attr("d",
               d3.svg.arc().innerRadius(function(d){return (d.index!=i)? r0:r0+10;}).outerRadius(
                function(d){return (d.index!=i)? r1:r1+10;}));


          groupText.transition().attr("x",function(d){  
            if(d.index==i){
              if ((d.startAngle+d.endAngle)/2 < 3.1415)
                return 20;
              else
                return -10;
            }else{
              return 8;
            }
          });

          popUpList = [];

          popUpList.push(parseInt(i));//for consistency between two modes

          popUpHierarchy();

          chordLayout.transition().attr("d",d3.svg.chord().radius(function(d){return (d.index!=i)? r0:r0+10;}));

          
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

            if(go_contents.length > 1){//only sort if two go and more exist
              sortGOcontent(go_contents);
            }
            

            $('#content').append(renderGOTerm(go_contents));
            $('.dropbtn').css("width",detailPanelWidth-50);

            $('.gene_name').click(function(){
                $('#filter').val($(this).html());
                refreshDetailPanel();
            });

            createOnClick(num_array);
            createOnHover(num_array);

            popUpList=[];
            for(i in index_array){
              var num = parseInt(index_array[i]);
              popUpList.push(num);
            }

            groupLayout.transition().attr("d",
                 d3.svg.arc().innerRadius(function(d){return (popUpList.indexOf(d.index)==-1)? r0:r0+5;}).outerRadius(
                  function(d,i){return (popUpList.indexOf(d.index)==-1)? r1:r1+5;}));

            groupText.transition().attr("x",function(d){
              if(popUpList.indexOf(d.index)!=-1){
                if ((d.startAngle+d.endAngle)/2 < 3.1415)
                  return 20;
                else
                  return -10;
              }else{
                return 8;
              }
            });

            popUpHierarchy();

            chordLayout.transition().attr("d",d3.svg.chord().radius(function(d){return (popUpList.indexOf(d.index)==-1)? r0:r0+5;}));

         /*   if(popUpList.length!=0)
              moveOutHierCluster();*/
        }
        else{
          resetVis();
          resetHierarchy();
          popUpList= [];

         
        }
      }
    }

    function popUpArcs(){
      clusterArc
        .style("fill", "green")
        .attr("class","clusterArc")
        .transition()
        .duration(300)
        .attr("d",function(d){
            
            var arc = d3.svg.arc()
              .innerRadius(d["radius"]+15)
              .outerRadius(d["radius"]+16)
              .startAngle(d["startAngle"])
              .endAngle(d["endAngle"]);

            return arc();
        });
    }

    function resetArcs(){
      clusterArc
        .style("fill", "green")
        .attr("class","clusterArc")
        .transition()
        .duration(300)
        .attr("d",function(d){
            
            var arc = d3.svg.arc()
              .innerRadius(d["radius"]+10)
              .outerRadius(d["radius"]+11)
              .startAngle(d["startAngle"])
              .endAngle(d["endAngle"]);

            return arc();
        });
    }

    function popUpLines(){
      clusterLine1
        .style("fill", "green")
        .attr("class","clusterLine")
        .transition()
        .duration(300)
        .attr("d",function(d){
            
            var firstLine = d3.svg.arc()
            .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+10:d["LineInnerPosition"])
            .outerRadius(d["LineOuterPosition"]+10)
            .startAngle(d["LineAngle"])
            .endAngle(d["LineAngle"]+0.002);

            return firstLine();
        });

      clusterLine2
        .style("fill", "green")
        .attr("class","clusterLine")
        .transition()
        .duration(300)
        .attr("d",function(d){

          var secondLine = d3.svg.arc()
            .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+10:d["LineInnerPosition"])
            .outerRadius(d["LineOuterPosition"]+10)
            .startAngle(d["LineAngle"])
            .endAngle(d["LineAngle"]+0.002);

            return secondLine();
        });
    }

    function resetLines(){
      clusterLine1
        .style("fill", "green")
        .attr("class","clusterLine")
        .transition()
        .duration(300)
        .attr("d",function(d){
            
            var firstLine = d3.svg.arc()
            .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+5:d["LineInnerPosition"])
            .outerRadius(d["LineOuterPosition"]+5)
            .startAngle(d["LineAngle"])
            .endAngle(d["LineAngle"]+0.002);

            return firstLine();
        });

      clusterLine2
        .style("fill", "green")
        .attr("class","clusterLine")
        .transition()
        .duration(300)
        .attr("d",function(d){

          var secondLine = d3.svg.arc()
            .innerRadius((d["index"]>nodesSize)?d["LineInnerPosition"]+5:d["LineInnerPosition"])
            .outerRadius(d["LineOuterPosition"]+5)
            .startAngle(d["LineAngle"])
            .endAngle(d["LineAngle"]+0.002);

            return secondLine();
        });
    }

    function popUpClusterNodes(){
      clusterHierNodeView.transition().duration(300).attrTween("transform", tween).attr("class","clusterNodeView");
      clusterHierNodeTextView.transition().duration(300).attrTween("transform", tween).attr("class","clusterText");

      function tween(d, i, a) {
              var interpolate;
              var str;
              
              str = "rotate(" + (d["angle"] * 180 / Math.PI - 90) + ")"
               + "translate(" + ( d["radius"] + 10 ) + ",0)";
              interpolate = d3.interpolate(a,str);


              return function(t) {
                  return interpolate(t);
              };
      }
    }

    function resetClusterNodes(){
      clusterHierNodeView.transition().duration(300).attrTween("transform", tween).attr("class","clusterNodeView");
      clusterHierNodeTextView.transition().duration(300).attrTween("transform", tween).attr("class","clusterText");

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

    function popUpHierarchy(){
      popUpArcs();
      popUpLines();
      popUpClusterNodes();
    }

    function resetHierarchy(){
      resetArcs();
      resetLines();
      resetClusterNodes();
    }

    function drawLable(){

      labelOff = 0;

      // clusterArc.style("display","none");
      // clusterLine1.style("display","none");
      // clusterLine2.style("display","none");
      // clusterHierNodeView.style("display","none");
      // clusterHierNodeTextView.style("display","none");
      // groupText.style("display","")

      updateLabel();
    }

    function updateLabel(){
      groupText.remove();

      groupText = circleSvg.append("svg:g")
            .selectAll("g")
            .data(chord.groups)
            .enter().append("svg:g")
            .attr("transform", function(d) {
            return "rotate(" + (((d.startAngle+d.endAngle)/2 * 180 / Math.PI - 90)+5) + ")"
               + "translate(" + (r1 + 0) + ",0)";
            }).append("svg:text")
            .style('fill', 'grey')
            .attr("x", function(d){
              if(popUpList.indexOf(d.index)!=-1){
                if ((d.startAngle+d.endAngle)/2 < 3.1415)
                  return 20;
                else
                  return -10;
              }else{
                return 8;
              }
            })
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

    function drawPopUpHierarchy(){
      setUpElementIfNecessary();
      popUpHierarchy();
    }

    function setUpElementIfNecessary(){
        drawHierarchicalClustering(clusterHierData);
    }

    function drawHierCluster(){
      labelOff = 1;
      // clusterArc.style("display","");
      // clusterLine1.style("display","");
      // clusterLine2.style("display","");
      // clusterHierNodeView.style("display","");
      // clusterHierNodeTextView.style("display","");


      // groupText.style("display","none");
      
      if(popUpList.length != 0){
        drawPopUpHierarchy();
      }else{
        drawHierarchicalClustering(that.clusterHierData);
      }
      

      if(textBackground){
        if(zoom.scale()<0.7){
            textBackground.attr("visibility","hidden");
        }else{
            textBackground.attr("visibility","visible");
        }
      }
    }

    function toggleDetails(){
      $('#arrow_detailedPanel').css('transform', function(){ return details_opened ? 'rotate(0deg)' : 'rotate(180deg)'});
      var shiftRight = detailPanelWidth - 50;

      $('#details').css('margin-right', function(){ return details_opened ? '-'+ shiftRight +'px' : '0'});
      details_opened = !details_opened;

      redrawMain_vis(details_opened);
      movePvalPanel(details_opened);
    }

    function toggleControl(){

      $('#arrow_controlPanel').css('transform', function(){ return control_opened ? 'rotate(180deg)' : 'rotate(0deg)'});
      
      var shiftleft = controlPanelWidth - 20;
      $('#control-panel').css('margin-left', function(){ return control_opened ? '-' + shiftleft+ 'px' : '0'});
      control_opened = !control_opened;
    }

    function redrawMain_vis(details_opened){

      if(!details_opened){
        d3.select(".main_vis").attr("width",w+detailPanelWidth);
        circleSvg.transition().attr("transform", "translate(" + (w+detailPanelWidth) / 2 + "," + h / 2 + ") scale(" + zoom.scale() + ")");
      }else{
        d3.select(".main_vis").attr("width",w);
        circleSvg.transition().attr("transform", "translate(" + w / 2 + "," + h / 2 + ") scale(" + zoom.scale() + ")");
      }
    }

    function movePvalPanel(details_opened){
      if(!details_opened){
        d3.select(".pval-label").transition().duration(300).style("margin-left",-210);
      }else{
        d3.select(".pval-label").transition().duration(300).style("margin-left",-detailPanelWidth-160);
      }
    }

    function setUpControlPanel(){

      if (reloadData == false){
          var element = '&nbsp<label>Cluster GO term according to the minimum number of common gene(s)</label> \
              <div id="slider" class="sliderBar"></div>\
              <input type="text" id="input_slider"/>';
                // <table class="RadioBox">\
                // <tbody><tr><td>\
                //       <input id="labelRadioBox" class="radioButton" type="radio" name="radioBox" value="0" checked>\
                //     </td><td><label style="padding-left:10px" for="labelRadioBox">Show GO term name</label>\
                //     </td></tr><tr><td></td></tr><tr><td>\
                //       <input id="hierClusterRadioBox"  class="radioButton" type="radio" name="radioBox" value="1" >\
                //     </td><td>\
                //        <label style="width:300px;padding-left:10px" for="hierClusterRadioBox">Show hierarchical tree and click on the node to manually cluster the GO term</label>\
                //     </td></tr></tbody></table>';

          //add save image button
          element += '<div class="control-panel-button">';
          element += '<button id="editor_save" class="btn" z-index:100">Save image</button>';
          element +=  '<button id="export" class="btn" z-index:100">Export File</button>';
          element +=  '<input type="file" id="import" style="display: none"> <label for="import" class="btn" id="import_label">Import File</label>'
          
          // element +=  '<input id="import" type="file" class="btn" z-index:100"></input>';
          element += '</div>';
          //add toggle button
          element +=  '<button id="arrow_controlPanel" class="arrowBar arrow"></button>';

        $("#control-panel").append(element);

      }

      var ranger = document.getElementById('slider');
      var inputFormat = document.getElementById('input_slider');
      setUpRangeSlider(ranger,inputFormat,that.minNumOfOverlappingGens,that.maxNumOfOverlappingGens);
    }

    function setUpRangeSlider(ranger,inputFormat,minNumOfOverlappingGens,maxNumOfOverlappingGens){
      if(ranger.noUiSlider)
        ranger.noUiSlider.destroy();

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
          var level = getLevelFromNumOfOverlappingGenes(values.toString());
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

    function getRidOfFadeLines(svg){
        return svg.replace(/<path class="chord fade" d="[,-\d\s\w\.]*" id="chordChords"><\/path>/g, "");
    }

    function setUpListener(monago){
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

      if (reloadData == false){
        $('#arrow_detailedPanel').click(function(){
          toggleDetails();
        });

        $('#arrow_controlPanel').click(function(){
          toggleControl();
        });  
      }


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

      $("#editor_save").click(function() {


          d3.select("canvas").attr("width",w).attr("height",h);

          var svg = d3.select('.main_vis')
                    .attr("version", 1.1)
                    .attr("xmlns", "http://www.w3.org/2000/svg")
                    .node().parentNode.innerHTML;


          svg = getRidOfFadeLines(svg);
          //console.log(svg);

          var imgsrc = 'data:image/svg+xml;base64,'+ btoa(svg);

          var canvas = document.querySelector("canvas");
          var context = canvas.getContext("2d");

          var image = new Image();
          image.src = imgsrc;

          image.onload = function() {
            //clean the context for redrawing

            context.clearRect(0, 0, w, h);

            context.drawImage(image, 0, 0);

            var canvasdata = canvas.toDataURL("image/png");
/*            var a = document.createElement("a");
            a.download = "sample.png";
            a.href = canvasdata;
            a.click();*/
            post("/getPic", {
                filename: 'chart',
                png: canvasdata
            }, null);

          };

        });


      $('#export').click(function(){

          var save_file = "{\"size\":" + size + "," + "\"go_inf\":" + JSON.stringify(go_inf) + "," + "\"clusterHierData\":" +JSON.stringify(clusterHierData)
          + "," + "\"matrix\":"+ JSON.stringify(that.matrix) + "," +"\"array_order\":" + JSON.stringify(array_order) + "," + "\"goNodes\":" + JSON.stringify(goNodes)+"}";

            post("/export", {
                filename: 'file',
                file: save_file
            }, null);
      });

      $('#import').click(function(){

            function readSingleFile(evt) {
            //Retrieve the first (and only!) File from the FileList object
            var f = evt.target.files[0]; 

            if (f) {
              var r = new FileReader();

              r.onload = function(e) { 
                var contents = e.target.result; 
                monago.reload(contents);
              }
              r.readAsText(f);
            } else { 
              alert("Failed to load file");
            }
          }

        document.getElementById('import').addEventListener('change', readSingleFile, false);
      })
    }

    function mextend (a, b) {
          var n;
          if (!a) {
              a = {};
          }
          for (n in b) {
              a[n] = b[n];
          }
          return a;
    };


    function mcss(el, styles) {
          // if (H.isMS && !H.svg) { // #2686
          //     if (styles && styles.opacity !== undefined) {
          //         styles.filter = 'alpha(opacity=' + (styles.opacity * 100) + ')';
          //     }
          // }
          mextend(el.style, styles);
      };

    function createElement(tag, attribs, styles, parent, nopad) {
            var el = window.document.createElement(tag),
                css = mcss;
            if (attribs) {
                mextend(el, attribs);
            }
            if (nopad) {
                css(el, {
                    padding: 0,
                    border: 'none',
                    margin: 0
                });
            }
            if (styles) {
                css(el, styles);
            }
            if (parent) {
                parent.appendChild(el);
            }
            return el;
    };

    // Add the H.post utility
    function post(url, data, formAttributes) {

        var name,
            form;

        // create the form
        form = createElement('form', {
            method: 'post',
            action: url,
            enctype: 'multipart/form-data'
        }, {
            display: 'none'
        }, window.document.body);

        // add the data
        for (name in data) {
            createElement('input', {
                type: 'hidden',
                name: name,
                value: data[name]
            }, null, form);
        }

        // submit
        form.submit();

        // clean up
        discardElement(form);
    };


    function discardElement(element) {
      element.remove();
    };


    function setUpView(){

      if(main_div)
        main_div.remove();

      main_div = d3.select("#chart").append("div");

      $("#details").css("width",detailPanelWidth);
      $("#pval-label").css("margin-left",-detailPanelWidth-160);
      $("#control-panel").css("width",controlPanelWidth);
      $("#control-panel").css("height",controlPanelHeight);

      svg = main_div
       .append("svg:svg")
         .attr("class","main_vis")
         .attr("width", w)
         .attr("height", h)
         .call(zoom.on("zoom", redraw))

      circleSvg = svg.append("svg:g")
         .attr("id", "circle")
         .attr("transform", "translate(" + w / 2 + "," + h / 2 + ")")

      circleSvg.append("circle").attr("r",r0).style("opacity", "0");
      
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

      groupText = circleSvg.append("svg:g")
           .selectAll("g")
             .data(chord.groups)
            .enter().append("svg:g");

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

    function copyGOInfFrom(go_inf){
      var go_inf_ori = [];
      go_inf.forEach(function(d,i){
        go_inf_ori.push(d);
      });
      return go_inf_ori;

    }

    function perparePvalue(){
      seperated_points = [];
      maxpVal = getMaxPval();
      minpVal = getMinPval();

      pVal_sort = sortGO_inf();

      step = (maxpVal-minpVal)/10;
      
      for(i=0;i<11;i++){
        seperated_points.push(minpVal+step*i);
      }
    }

    function perpareChord(){
      chordMatrix = d3.layout.chord()
       .padding(.03)
       .sortSubgroups(d3.descending);
       

      chord = chordMatrix.matrix(that.matrix);

      //store the nodes for hierarchical clustering visualzaition
      chordGroupsNodePosition = [];
      chord.groups().forEach(function(d,i){
        nodeObj = {"index":d.index,"angle":(d.startAngle+d.endAngle)/2,"radius":r1};
        chordGroupsNodePosition.push(nodeObj);
      });
    }

    this.init = function(size,go_inf,clusterHierData,goNodes,matrix){

      that.go_inf = go_inf;
      that.goNodes = goNodes;
      that.go_inf_ori = copyGOInfFrom(go_inf);
      that.matrix = matrix;

      that.clusterHierData = clusterHierData;


      that.maxNumOfOverlappingGens = that.clusterHierData[0][3];
      that.minNumOfOverlappingGens = that.clusterHierData[that.clusterHierData.length-1][3];

      that.clusterHierData.map(function(d){
        that.clusterHierDataStatic.push([d[2],d[3]]);
      });

      //get a list of genes for each GO term
      that.go_inf.forEach(function(d,i){
        that.go_inf[i].genes = d.genes.split("|");
      });

      perparePvalue();

      perpareChord();

      createGeneInfo();
      createInteractionGenes();

      createPvalLabelPanel(fill);
      createClusterNodesStatus(clusterHierData,[],"",[]);
      setUpView();
      determineLabelSize();
      drawHierCluster();
      drawLable();

      setUpListener(that);
      
    }

    this.reload = function(content){
      reloadData = true;
      that.level_g = 0;

      var content = JSON.parse(content);

      that.go_inf = content["go_inf"];
      that.goNodes = content["goNodes"];
      that.go_inf_ori = copyGOInfFrom(that.go_inf);
      that.matrix = content["matrix"];


      array_order = content["array_order"];

      nodesSize = content["size"];
      maxNodeIndex = content["size"];

      that.clusterHierData = content["clusterHierData"];
      that.maxNumOfOverlappingGens = that.clusterHierData[0][3];
      that.minNumOfOverlappingGens = that.clusterHierData[that.clusterHierData.length-1][3];


      that.clusterHierDataStatic = [];
      clusterHierNodesStatus = [];
      that.memCache = {};

      that.clusterHierData.map(function(d){
        that.clusterHierDataStatic.push([d[2],d[3]]);
      });

      perparePvalue();

      perpareChord();

      createGeneInfo();
      createInteractionGenes();

      createPvalLabelPanel(fill);
      createClusterNodesStatus(that.clusterHierData,[],"",[]);
      setUpView();
      determineLabelSize();

      clusterArc.remove()
      clusterLine2.remove()
      clusterLine1.remove()

      drawHierCluster();
      drawLable();

      setUpListener();

    }

  }

  monago = new MonaGO().init(size,go_inf,clusterHierData,goNodes,matrix);


})();










