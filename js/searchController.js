var app = angular.module('MonaGO', []);

app.constant('goNodes', goNodes);

app.controller('searchCtrl', ['$scope', 'goNodes', function($scope,goNodes) {
    $scope.goNodes = [];

    for (var property in goNodes) {
	    if (goNodes.hasOwnProperty(property)) {

	        $scope.goNodes.push({"id":property,"term":goNodes[property]['n']});
	    }
	}	
	console.log(goNodes);
    
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
    	
    }
}]);