var app = angular.module('MonaGO', []);

app.constant('go_inf', go_inf);

app.controller('searchCtrl', ['$scope', 'go_inf', function($scope,go_inf) {
    $scope.goArr = [];

/*    for (var property in go_inf) {
	    if (go_inf.hasOwnProperty(property)) {

	        $scope.go_inf.push({"GO_name":property,"GO_id":goNodes[property]['GO_id']});
	    }
	}	*/

    go_inf.forEach(function(d,i){
        $scope.goArr.push({"id":d.GO_id,"name":d.GO_name,"genes":d.genes});
    });
    
    $scope.isNotEmpty = function(searchText){

        if (searchText.length!=0){
            return true;
        }
        else{
            return false;
        }
    }

    $scope.setSearchBox = function(GO_id){
    	$scope.searchText = GO_id;
        $('#searchBox').remove();
    	
    }
}]);