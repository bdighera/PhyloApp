// /**
//  * Created by Bryan on 2/20/21.
//  */
function getRequest(input){

	var XMLHttpRequest = require("xmlhttprequest").XMLHttpRequest;

// initializes and establishes url for get request
const xhr = new XMLHttpRequest();
var base = 'http://127.0.0.1:5000/domains/'
var sendrequest = base.concat(input);
const url = sendrequest;

//sets response type and parses return if successful request
xhr.responseType = 'json';
xhr.onreadystatechange = function() {
  if (this.status == 200 && xhr.readyState == 4)
      var json = JSON.parse(this.responseText);
        array = RemoveUndefined(json);


  };

// removes the undefined elements in the json
function RemoveUndefined(json) {
    if (json != undefined){return console.log(json)}}

xhr.open('GET', url);
xhr.send();
}


function parser(){}

validate = function() {
    message = document.getElementById("name").value;
    correctMessage = message.replace('%2C',',')
    //console.log(correctMessage)

}

function getWholeDB(){

    	var XMLHttpRequest = require("xmlhttprequest").XMLHttpRequest;

        // initializes and establishes url for get request
        const xhr = new XMLHttpRequest();
        var base = 'http://127.0.0.1:5000//DataBaseDisplay/';
        // var sendrequest = base.concat(input);
        //const url = sendrequest;

        //sets response type and parses return if successful request
        xhr.responseType = 'json';
        xhr.onreadystatechange = function() {
          if (this.status == 200 && xhr.readyState == 4)
              var json = JSON.parse(this.responseText);
                array = RemoveUndefined(json);


  };
}

getWholeDB();

//NP_001243399.1,XP_016778522.1,XP_016874984.1


