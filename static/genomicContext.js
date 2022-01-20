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
var domainCount = 0;
var legend = [];
var geneLegend = [];
var domains = [];
var genes = [];
var allDomains = [];
var allMotifs = [];
var selectedMotifs = [];
var colorStack = ['yellow', 'green'];
var highlightedMotifs = 0;

var newData;

function postSequences() {
    let postData = '';
    for(let sequence of deletedSequences){
        postData += sequence + ' ';
    }
    alert(deletedSequences);
    $.post("http://localhost:8080/genomicContext", {'deleted_sequences':postData});
}

function compareMotifs() {
  if(selectedMotifs.length == 2){
    let postData = [];
    for(let motif of selectedMotifs){
      postData.push(motif);
    }
    alert(postData);
    $.post("http://localhost:8080/GCAlign", {'compared_motifs':postData});
  }else{
    alert("Please highlight two motifs");
  }
}

function sequenceData(data) {
	let sequences = data.Sequences;
	var sequenceId = 0;
  newData = data
	for(let sequence of sequences) {
    var newSequence = new Sequence(sequence, sequenceId);
		sequenceId += 1;
	}
	displayLegend();
}

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
    this.parseIndex = data.name.indexOf('.');
    if (this.parseIndex != -1) {
		  this.name = data.name.substring(0, this.parseIndex);
    } else {
		  this.name = data.name;
    }
		this.sequenceId = sequenceId;
		this.motifId = motifId;
		this.domainId = domainId;
    this.color = data.color;
	}

	addColor(){
		let fullId = [this.sequenceId, this.motifId, this.domainId];
		if(domains.includes(this.name, 0)) {
			var thing = legend.find( ({name}) => name === this.name);
			thing.addMotifId(fullId);
      this.color = thing.color;
      newData.Sequences[this.sequenceId].motifs[this.motifId].domains[this.domainId].color = thing.color;
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

    var genomicContextTable = document.getElementById('sequenceTable');

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
    let sequenceNameText = document.createElement("h3");
	  sequenceNameText.innerHTML += this.sequenceName;
    sequenceNameText.style.marginBottom = "0px";
    sequenceNameDiv.appendChild(sequenceNameText);

    let optionsDiv = document.createElement("div");
	  let reverseDiv = document.createElement("div");
    let reverseInput = document.createElement("input");
	  reverseInput.type = "checkbox";
	  reverseInput.addEventListener('change', onClickReverse);
	  reverseInput.name = this.sequenceName;
    let reverseText = document.createElement("label");
    reverseText.textContent = "Reverse";
    reverseDiv.appendChild(reverseInput);
    reverseDiv.appendChild(reverseText);

    let deleteDiv = document.createElement("div");
    let deleteInput = document.createElement("input");
    deleteInput.type = "checkbox";
    deleteInput.name = this.sequenceName;
    deleteInput.addEventListener('change', onClickDelete);
    let deleteText = document.createElement("label");
    deleteText.textContent = "Delete";
    deleteDiv.appendChild(deleteInput);
    deleteDiv.appendChild(deleteText);
    deleteDiv.style.width="50%";
    reverseDiv.style.width="50%";
    optionsDiv.appendChild(reverseDiv);
    optionsDiv.appendChild(deleteDiv);
    optionsDiv.style.display="flex";

    this.optionsCell.appendChild(sequenceNameDiv);
    this.optionsCell.appendChild(optionsDiv);
  }

  addSequence(){
    this.sequenceCell = this.sequenceRow.insertCell(); 
	  var motifId = 0;
	  for(let motif of this.motifs) {
      let newMotif = new MotifDivNode(this.sequenceDiv, motifId, motif, this.id, this.sequenceName);
      allMotifs.push(newMotif);
		  motifId += 1;
	  }

    let maxOrder = motifId*2;
    let startDiv = new EmptyDivNode(this.sequenceDiv, 0, maxOrder);
    let endDiv = new EmptyDivNode(this.sequenceDiv, maxOrder, maxOrder);
  }
}

function addDomain(ctx, domain, motifName, sequenceId, motifId, domainId, domainLength) {
	let start = domainId * domainLength;

	var newDomain = new Domain(domain, sequenceId, motifId, domainId);
	newDomain.addColor();
	allDomains.push(newDomain);
	domainCount += 1;
  
	ctx.fillStyle = newDomain.color;
  ctx.fillRect(start + domainStart, domainTop, domainLength, domainHeight);
  return newDomain;
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

function editGene(motifId, motif, chosenColor) {
	var canvas = document.getElementById(motifId);
	var ctx = canvas.getContext('2d');
	ctx.clearRect(0,0, motifWidth, motifHeight);

	domainBox(ctx);
	direction(ctx, motif.direction);	

  ctx.font = "24px Helvetica";
  ctx.fillStyle = chosenColor;
  ctx.fillRect(domainStart, domainTop, domainEnd-arrowWidth-domainTop, domainHeight);
  ctx.fillStyle = "#000";
  ctx.fillText(motif.geneName, 26, 31, domainWidth);
}

function editMotif(motifId, motif) {
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

class DivNode {
	constructor(order, sequence) {
    this.div = document.createElement("div");
    this.div.style.width = "250";
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
class MotifColor {
	constructor(geneName, geneColor, motifId) {
		this.name = geneName;
		this.color = geneColor;
		this.motifIds = [motifId];
	}

	addToGeneLegend() {
		genes.push(this.name);
		geneLegend.push(this);
	}

	addMotifId(motifId) {
		this.motifIds.push(motifId);
	}
}

class MotifDivNode extends DivNode {
	constructor(sequence, order, motif, sequenceId, sequenceName) {
		super(sequence, order);
	
    this.motifId = order;
    this.order = 1+(2*order);
		this.motif = motif;
		this.canvas = document.createElement('canvas');
		this.context = this.canvas.getContext("2d");
    this.sequenceName = sequenceName;
    this.geneName = this.motif.geneName;
    this.sequence = this.motif.seq

    this.motifName = this.sequenceName + this.geneName;
    const motifName = [this.sequenceName,this.geneName];

    this.canvas.height = motifHeight;
    this.canvas.width = motifWidth;	
    this.canvas.id = this.motifName;

  
    this.fullId = [sequenceId, this.motifId];
    if(genes.includes(this.geneName, 0)) {
      var thing = geneLegend.find( ({name}) => name === this.geneName.toUpperCase());
      thing.addMotifId(this.fullId);
    } else {
      let newMotif = new MotifColor(this.geneName.toUpperCase(), this.color, this.fullId);
      newMotif.addToGeneLegend();
    }

    var _this = this;

    this.div.addEventListener('click', function (event) {
      if(this.style.backgroundColor != 'yellow' && this.style.backgroundColor != 'green' && highlightedMotifs < 2) {
        var newColor = colorStack.pop();
        for(let motif of allMotifs) {
          if(motif.geneName.toUpperCase() == _this.geneName.toUpperCase()) {
              motif.div.style.backgroundColor = newColor;
            
          }
        }
        highlightedMotifs += 1;
        selectedMotifs.push(_this.sequence.toUpperCase());
      } else if(this.style.backgroundColor == 'yellow' || this.style.backgroundColor == 'green') {
        const removedMotifName = selectedMotifs.indexOf(_this.sequence.toUpperCase());
        colorStack.push(this.style.backgroundColor);
        if(removedMotifName > -1) {
          selectedMotifs.splice(removedMotifName,1);
        }
        for(let motif of allMotifs) {
          if(motif.geneName.toUpperCase() == _this.geneName.toUpperCase()) {
            motif.div.style.backgroundColor = 'transparent';
          }
        }
        highlightedMotifs -= 1;
      }
      
    });

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
        let newDomain = addDomain(this.context, domain, this.motifName, sequenceId, this.motifId, domainId, this.domainLength);
        domainId += 1;
        domainList += '\t' + newDomain.name + '\n';
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

class EmptyDivNode extends DivNode {
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
      let newNode = new EmptyDivNode(sequence, parseInt(this.div.style.order), this.maxOrder)
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

function displayLegend() {
  var legendTable = document.getElementById('geneLegend');
  for(let item of geneLegend) {
		let legendRow = legendTable.insertRow();

		let legendNameCell = legendRow.insertCell();
		let legendNameDiv = document.createElement("div");
		legendNameDiv.style.textAlign="right";
		legendNameDiv.innerHTML += item.name;	
		legendNameCell.appendChild(legendNameDiv);
    
		let legendColorCell = legendRow.insertCell();
		let changeColorDiv = document.createElement("input");
		changeColorDiv.type = "color";
    changeColorDiv.value = "#FFFFFF";
		changeColorDiv.name = item.name;
		changeColorDiv.addEventListener('change', function() {
      changeGeneColor(changeColorDiv.value, item.motifIds)
    });
		legendColorCell.appendChild(changeColorDiv);
  }

  var legendTable = document.getElementById('legend');
	for(let item of legend) {
		let legendRow = legendTable.insertRow();

		let legendColorCell = legendRow.insertCell();
		let changeColorDiv = document.createElement("input");
		changeColorDiv.type = "color";
    changeColorDiv.value = item.color;
		changeColorDiv.name = item.name;
		changeColorDiv.addEventListener('change', function() {
      changeColor(changeColorDiv.value, item.motifIds)
    });
		legendColorCell.appendChild(changeColorDiv);

		let legendNameCell = legendRow.insertCell();
		let legendNameDiv = document.createElement("div");
		legendNameDiv.style.textAlign="left";
		legendNameDiv.innerHTML += item.name;	
		legendNameCell.appendChild(legendNameDiv);

	}
}

function changeGeneColor(newColor, domains) {

	for(let domain of domains) {
		let motif = newData.Sequences[domain[0]].motifs[domain[1]];
		let motifId = newData.Sequences[domain[0]].name + newData.Sequences[domain[0]].motifs[domain[1]].geneName;
		newData.Sequences[domain[0]].motifs[domain[1]].color = newColor;
		editGene(motifId, motif, newColor);
  }
}
  

function changeColor(newColor, domains) {
	for(let domain of domains) {
		let motif = newData.Sequences[domain[0]].motifs[domain[1]];
		let motifId = newData.Sequences[domain[0]].name + newData.Sequences[domain[0]].motifs[domain[1]].geneName;
		newData.Sequences[domain[0]].motifs[domain[1]].domains[domain[2]].color = newColor;
		editMotif(motifId, motif);
  }
}

