var app = angular.module('MonaGO', []);

app.constant('goNodes', goNodes);

app.controller('searchCtrl', ['$scope', 'goNodes', function($scope,goNodes) {
    $scope.goNodes = goNodes;

    $scope.isNotEmpty = function(searchText){

        if (searchText.length!=0){
            return true;
        }
        else{
            return false;
        }
    }
}]);