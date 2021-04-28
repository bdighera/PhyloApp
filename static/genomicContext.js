$(document).ready(function(){


var deletedSequences = [];

var domainHeight = 48;
var domainTop = 1;
var domainStart = 26;
var domainEnd = 225;
var domainWidth = domainEnd - domainStart;
var motifEnd = 250;
var motifStart = 0;
var motifTop = 0;
var motifBottom = 50;
var motifHeight = 50;
var motifWidth = 250;
var arrowWidth = 25;

var allNodes = [];

var genomicContextTable = document.getElementById('sequenceTable');
var genomicContextDiv = document.getElementById('demo');
var legendTable = document.getElementById('legend');

var exportButton = document.getElementById('exportButton');
//TODO: Need to get the export button working.
var domainCount = 0;
var legend = [];
var domains = [];
var allDomains = [];

var motifList = [];

var lastDomain = [];
// var stuff;
// $.get( "http://localhost:5000/genomicContext.html", function( gcData ) {
// 	$( ".result" ).html( data );
// 	stuff = data;
// 	sequenceData(data);
// });

class DomainColor {
	constructor(domainName, domainColor, motifId) {
		this.name = domainName;
		this.color = domainColor;
		this.motifIds = [motifId];
	}

	addToLegend() {
		domains.push(this.name);
		legend.push(this);
	}

	addMotifId(motifId) {
		this.motifIds.push(motifId);
	}
}

class Domain {
	constructor(data, sequenceId, motifId, domainId) {
		this.name = data.name;
		this.color = data.color;
		this.sequenceId = sequenceId;
		this.motifId = motifId;
		this.domainId = domainId;
	}

	addColor(){
		lastDomain = this.domainId;
		let fullId = [this.sequenceId, this.motifId, this.domainId];
		if(domains.includes(this.name, 0)) {
			var thing = legend.find( ({name}) => name === this.name);
			thing.addMotifId(fullId);
		} else {
			let newDomain = new DomainColor(this.name, this.color, fullId);
			newDomain.addToLegend();
		}
	}
}

class Sequence {
	constructor(data, id) {
    this.id = id;
		this.sequenceName = data.name;
		this.motifs = data.motifs;
    this.divNodes = [];
    this.divNodesOrder = "forward";

    this.sequenceRow = genomicContextTable.insertRow();
    
    this.sequenceDiv = document.createElement("div");
    this.sequenceDiv.id = this.sequenceName + "div";
    this.sequenceDiv.className = "sequence";

    this.addOptions();
  
    this.addSequence();  

	  this.sequenceCell.appendChild(this.sequenceDiv);	
	}

  addOptions(){
    this.optionsCell = this.sequenceRow.insertCell();
	  
    let sequenceNameDiv = document.createElement("div");
	  sequenceNameDiv.innerHTML += this.sequenceName;	

	  let reverseDiv = document.createElement("input");
	  reverseDiv.type = "checkbox";
	  reverseDiv.name = this.sequenceName;
	  reverseDiv.addEventListener('change', onClickReverse);

    let deleteDiv = document.createElement("input");
    deleteDiv.type = "checkbox";
    deleteDiv.name = this.sequenceName;
    deleteDiv.addEventListener('change', onClickDelete);

    this.optionsCell.appendChild(sequenceNameDiv);
    this.optionsCell.appendChild(reverseDiv);
    this.optionsCell.appendChild(deleteDiv);
  }

  addSequence(){
    this.sequenceCell = this.sequenceRow.insertCell(); 
	  var motifId = 0;
	  for(let motif of this.motifs) {
      let newMotif = new motifDivNode(this.sequenceDiv, motifId, motif, this.id, this.sequenceName);
		  motifId += 1;
	  }

 	  var motifOrder = this.divNodes;
    let maxOrder = motifId*2;
    let startDiv = new emptyDivNode(this.sequenceDiv, 0, maxOrder);
    let endDiv = new emptyDivNode(this.sequenceDiv, maxOrder, maxOrder);
  }
}

function addDomain(ctx, domain, motifName, sequenceId, motifId, domainId, domainLength) {
	let start = domainId * domainLength;

	var newDomain = new Domain(domain, sequenceId, motifId, domainId);
	newDomain.addColor();
	allDomains.push(newDomain);
	domainCount += 1;
  
	ctx.fillStyle = domain.color;
  ctx.fillRect(start + domainStart, domainTop, domainLength, domainHeight);

  // Domain name text
  //ctx.font = "9px Helvetica";
  //ctx.fillStyle = "#000";
  //ctx.fillText(domain.name, start + domainStart, 40, domainLength);
}

function direction(ctx, direction) {
  if (direction == ">") {
    ctx.fillStyle = "black";
    ctx.moveTo(motifEnd, motifHeight / 2);
    ctx.lineTo(motifEnd - arrowWidth, motifTop);
    ctx.lineTo(motifEnd - arrowWidth, motifBottom);
    ctx.fill();
  } else if (direction == "<") {
    ctx.fillStyle = "black";
    ctx.moveTo(motifStart, motifHeight / 2);
    ctx.lineTo(motifStart + arrowWidth, motifTop);
    ctx.lineTo(motifStart + arrowWidth, motifBottom);
    ctx.fill();
  }
}

function domainBox(ctx) {
  ctx.strokeRect(domainStart, motifTop, domainWidth, motifBottom);
}

function editMotif(motifId, motif, ids) {
	var canvas = document.getElementById(motifId);
	var ctx = canvas.getContext('2d');
	ctx.clearRect(0,0, motifWidth, motifHeight);

	let motifName = motifId;
	domainBox(ctx);
	direction(ctx, motif.direction);	

  var domainLength = domainWidth/motif.domains.length	
	var domainId = 0;
	for (let domain of motif.domains) {
		if (domain.name != null) {
			let start = domainId * domainLength;
			let end = 20;
  
			ctx.fillStyle = domain.color;
  		ctx.fillRect(start + domainStart, domainTop, domainLength, domainHeight);

      // Domain Name Text
  		//ctx.font = "9px Helvetica";
  		//ctx.fillStyle = "#000";
  		//ctx.fillText(domain.name, start + domainStart, 40, domainLength);
			domainId += 1;
		}
	}
  ctx.font = "24px Helvetica";
  ctx.fillStyle = "#000";
  ctx.fillText(motif.geneName, 26, 31, domainWidth);
}

function onClickReverse() {
	let sequenceId = this.name + "div"
	let sequenceDiv = document.getElementById(sequenceId);

	if (this.checked) {
		sequenceDiv.className = "reverseSequence";
	} else {
		sequenceDiv.className = "sequence";
	}
}

function onClickDelete() {
	let sequenceId = this.name + "div"
	let sequenceDiv = document.getElementById(sequenceId);

	if (this.checked) {
		sequenceDiv.style.opacity = "0.25";
		deletedSequences.push(this.name);
	} else {
		sequenceDiv.style.opacity = "1";
		let deletedLocation = deletedSequences.indexOf(this.name);
		if (deletedLocation > -1) {
			deletedSequences.splice(deletedLocation, 1);
		}
	}
}

class divNode {
	constructor(order, sequence) {
    this.div = document.createElement("div");
    this.div.style.width = "250";//motifEnd;
    this.order = parseInt(order);
		this.div.style.order = parseInt(order);
		this.sequence = sequence;
	}

	getOrder() {
		return parseInt(this.div.style.order);
	}

	incrementOrder() {
		this.div.style.order = this.getOrder() + 2;
	}

	decrementOrder() {
		this.div.style.order = this.getOrder() - 2;
	}
}

class motifDivNode extends divNode {
	constructor(sequence, order, motif, sequenceId, sequenceName) {
		super(sequence, order);
	
    this.motifId = order;
    this.order = 1+(2*order);
		this.motif = motif;
		this.canvas = document.createElement('canvas');
		this.context = this.canvas.getContext("2d");
    this.sequenceName = sequenceName;

    this.motifName = this.sequenceName + this.motif.geneName;

    this.canvas.height = motifHeight;
    this.canvas.width = motifWidth;	
    this.canvas.id = this.motifName;

    this.div.id = this.motifName + 'div';
    this.div.style.order = this.order;
    this.div.appendChild(this.canvas);	

    domainBox(this.context);
    direction(this.context, this.motif.direction);	
    this.domainLength = domainWidth/(this.motif.domains.length);
    var domainId = 0;
    var domainList = this.motif.geneName + " domains:\n";
    for (let domain of this.motif.domains) {
      if (domain.name != null) {
        addDomain(this.context, domain, this.motifName, sequenceId, this.motifId, domainId, this.domainLength);
        domainId += 1;
        domainList += '\t' + domain.name + '\n';
      }
    }

    if (domainId == 0) { 
      domainList += '\tnone';
    }
   
    this.div.title = domainList;
    this.context.font = "24px Helvetica";
    this.context.fillStyle = "#000";
    this.context.fillText(motif.geneName, 26, 31, domainWidth);
    this.div.style.width="500";
    sequence.appendChild(this.div);
	}
}

class emptyDivNode extends divNode {
	constructor(sequence, order, maxOrder) {
		super(sequence, order);   
    this.sequence = sequence;
    this.order = order;
    this.maxOrder = maxOrder;
    this.div.style.order=order;
    this.div.style.width="250px";
    this.div.style.display="flex";
    this.div.style.alignItems="center";
    this.div.style.justifyContent="center";
    var motifDivs = sequence.childNodes;

    this.leftButton = document.createElement("input");
    this.leftButton.type = "button";
    this.leftButton.value = "←";
    this.leftButton.addEventListener('click', ()=>{
      this.moveLeft();
    });
    this.div.appendChild(this.leftButton);

    this.delDivButton = document.createElement("input");
    this.delDivButton.type = "button";
    this.delDivButton.value = "-";
    this.delDivButton.addEventListener('click', ()=>{
      this.div.remove();
    });
    this.div.appendChild(this.delDivButton);

    this.addDivButton = document.createElement("input");
    this.addDivButton.type = "button";
    this.addDivButton.value = "+";
    this.addDivButton.addEventListener('click', ()=>{
      let newNode = new emptyDivNode(sequence, parseInt(this.div.style.order), this.maxOrder)
    });
    this.div.appendChild(this.addDivButton);

    this.rightButton = document.createElement("input");
    this.rightButton.type = "button";
    this.rightButton.value = "→";
    this.rightButton.addEventListener('click', ()=>{
      this.moveRight();
    });
    this.div.appendChild(this.rightButton);

    sequence.appendChild(this.div);
	}

  addNewThing() {  
      let newNode = new emptyDivNode(sequence, this.order)
  }

	moveLeft(){
    let otherNode = this.sequence.childNodes[this.order];
    
    if(this.sequence.className == "sequence" && this.div.style.order > 0){
		  this.decrementOrder();
    } else if (this.sequence.className == "reverseSequence" && this.div.style.order < this.maxOrder){
      this.incrementOrder();
    }
	}	
	
	moveRight() {	
    let otherNode = this.sequence.childNodes[this.order];
    if(this.sequence.className == "sequence" && this.div.style.order < this.maxOrder) {
		  this.incrementOrder();
    } else if (this.sequence.className == "reverseSequence" && this.div.style.order > 0) {
      this.decrementOrder();
    }
	}
}

function DisplayFunct(data) {
	let sequences = data.Sequences;
	var sequenceId = 0;
	for(let sequence of sequences) {
    var newSequence = new Sequence(sequence, sequenceId);
		sequenceId += 1;
	}
	displayLegend();
}

function displayLegend() {
	for(let item of legend) {
		let legendRow = legendTable.insertRow();

		let legendNameCell = legendRow.insertCell();
		let legendNameDiv = document.createElement("div");
		legendNameDiv.style.textAlign="right";
		legendNameDiv.innerHTML += item.name;	
		legendNameCell.appendChild(legendNameDiv);

    let legendSwatchCell = legendRow.insertCell();
    let legendSwatchDiv = document.createElement("div");
    var  legendSwatchCanvas = document.createElement("canvas");
    legendSwatchCanvas.height = "25";
    legendSwatchCanvas.width = "50";	
    legendSwatchCanvas.id = item.color;
		var legendSwatchCtx = legendSwatchCanvas.getContext("2d");
    legendSwatchCtx.fillStyle = item.color;
    legendSwatchCtx.fillRect(0,0,100,50);
    legendSwatchDiv.appendChild(legendSwatchCanvas);
    legendSwatchCell.appendChild(legendSwatchDiv);	

		let legendColorCell = legendRow.insertCell();
		let changeColorDiv = document.createElement("input");
		changeColorDiv.type = "text";
		changeColorDiv.name = item.name;
		//changeColorDiv.value = item.color;
		changeColorDiv.addEventListener('change', function() {
      var newLegendSwatchCanvas = document.getElementById(item.color);
      var newLegendSwatchCtx = newLegendSwatchCanvas.getContext("2d");
      newLegendSwatchCtx.fillStyle = changeColorDiv.value;
      newLegendSwatchCtx.fillRect(0,0,100,50);
      changeColor(changeColorDiv.value, item.motifIds)
    });
		legendColorCell.appendChild(changeColorDiv);
	}
}

function changeColor(newColor, domains) {
	for(let domain of domains) {
		let motif = stuff.Sequences[domain[0]].motifs[domain[1]];
		let motifId = stuff.Sequences[domain[0]].name + stuff.Sequences[domain[0]].motifs[domain[1]].geneName;
		stuff.Sequences[domain[0]].motifs[domain[1]].domains[domain[2]].color = newColor;
		editMotif(motifId, motif, domain);
  }
}

function postSequences() {
	alert(deletedSequences);
}
});
