$(document).ready(function(){ 

var domainHeight = 48;
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

var genomicContextTable = document.getElementById('sequenceTable');
var genomicContextDiv = document.getElementById('demo');

$.get( "http://localhost:5000/genomicContext", function( data ) {
	$( ".result" ).html( data );

	sequenceData(data);		
});

function addDomain(ctx, domain, position) {
	let start = position * 20;
	let end = 20;
  ctx.fillStyle = domain.color;
  ctx.fillRect(start + domainStart, domainTop, end, domainHeight);

  ctx.font = "9px Helvetica";
  ctx.fillStyle = "#000";
  ctx.fillText(domain.name, start + domainStart, 40, end);
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

function addMotif(sequence, motif) {
	let motifDiv = document.createElement("div");
  let motifCanvas = document.createElement('canvas');
  motifCanvas.height = motifHeight;
  motifCanvas.width = motifWidth;
  motifCanvas.id = motif.name + 'Canvas';
	motifDiv.appendChild(motifCanvas);	

  ctx = motifCanvas.getContext("2d");
  domainBox(ctx);
  direction(ctx, motif.direction);

	var position = 0;
	for (let domain of motif.domains) {
		if (domain.name != null) {
  		addDomain(ctx, domain, position);
			position += 1;
		}
	}

  ctx.font = "24px Helvetica";
  ctx.fillStyle = "#000";
  ctx.fillText(motif.geneName, 26, 31, domainWidth);
	sequence.appendChild(motifDiv);
}

function addOptions(sequenceRow, sequence) {
	let optionsCell = sequenceRow.insertCell();
	let sequenceNameDiv = document.createElement("div");
	sequenceNameDiv.innerHTML += sequence.name;	

	let reverseDiv = document.createElement("input");
	reverseDiv.type = "checkbox";
	reverseDiv.name = sequence.name;
	reverseDiv.addEventListener('change', onClickReverse);

	let deleteDiv = document.createElement("input");
	deleteDiv.type = "checkbox";
	deleteDiv.name = sequence.name;
	deleteDiv.addEventListener('change', onClickDelete);

	optionsCell.appendChild(sequenceNameDiv);
	optionsCell.appendChild(reverseDiv);
	optionsCell.appendChild(deleteDiv);
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
	} else {
		sequenceDiv.style.opacity = "1";
	}
}

function addSequence(genomicContext, sequence) {
  let sequenceRow =  genomicContext.insertRow();
	let sequenceDiv = document.createElement("div");
	sequenceDiv.id = sequence.name + "div";
	addOptions(sequenceRow, sequence);
	sequenceDiv.className = "sequence";
	let motifs = sequence.motifs;
	let sequenceCell = sequenceRow.insertCell();
	
	for(let motif of motifs) {
		addMotif(sequenceDiv, motif);
	}
	sequenceCell.appendChild(sequenceDiv);	
}


function sequenceData(data) {
	let sequences = data.Sequences;
	for(let sequence of sequences) {
		let sequenceName = sequence.name;
		let motifs = sequence.motifs;
		addSequence(genomicContextTable, sequence);		
	}
}

});$(document).ready(function(){ 

var domainHeight = 48;
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

var genomicContextTable = document.getElementById('sequenceTable');
var genomicContextDiv = document.getElementById('demo');

$.get( "http://localhost:5000/genomicContext", function( data ) {
	$( ".result" ).html( data );

	sequenceData(data);		
});

function addDomain(ctx, domain, position) {
	let start = position * 20;
	let end = 20;
  ctx.fillStyle = domain.color;
  ctx.fillRect(start + domainStart, domainTop, end, domainHeight);

  ctx.font = "9px Helvetica";
  ctx.fillStyle = "#000";
  ctx.fillText(domain.name, start + domainStart, 40, end);
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

function addMotif(sequence, motif) {
	let motifDiv = document.createElement("div");
  let motifCanvas = document.createElement('canvas');
  motifCanvas.height = motifHeight;
  motifCanvas.width = motifWidth;
  motifCanvas.id = motif.name + 'Canvas';
	motifDiv.appendChild(motifCanvas);	

  ctx = motifCanvas.getContext("2d");
  domainBox(ctx);
  direction(ctx, motif.direction);

	var position = 0;
	for (let domain of motif.domains) {
		if (domain.name != null) {
  		addDomain(ctx, domain, position);
			position += 1;
		}
	}

  ctx.font = "24px Helvetica";
  ctx.fillStyle = "#000";
  ctx.fillText(motif.geneName, 26, 31, domainWidth);
	sequence.appendChild(motifDiv);
}

function addOptions(sequenceRow, sequence) {
	let optionsCell = sequenceRow.insertCell();
	let sequenceNameDiv = document.createElement("div");
	sequenceNameDiv.innerHTML += sequence.name;	

	let reverseDiv = document.createElement("input");
	reverseDiv.type = "checkbox";
	reverseDiv.name = sequence.name;
	reverseDiv.addEventListener('change', onClickReverse);

	let deleteDiv = document.createElement("input");
	deleteDiv.type = "checkbox";
	deleteDiv.name = sequence.name;
	deleteDiv.addEventListener('change', onClickDelete);

	optionsCell.appendChild(sequenceNameDiv);
	optionsCell.appendChild(reverseDiv);
	optionsCell.appendChild(deleteDiv);
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
	} else {
		sequenceDiv.style.opacity = "1";
	}
}

function addSequence(genomicContext, sequence) {
  let sequenceRow =  genomicContext.insertRow();
	let sequenceDiv = document.createElement("div");
	sequenceDiv.id = sequence.name + "div";
	addOptions(sequenceRow, sequence);
	sequenceDiv.className = "sequence";
	let motifs = sequence.motifs;
	let sequenceCell = sequenceRow.insertCell();
	
	for(let motif of motifs) {
		addMotif(sequenceDiv, motif);
	}
	sequenceCell.appendChild(sequenceDiv);	
}


function sequenceData(data) {
	let sequences = data.Sequences;
	for(let sequence of sequences) {
		let sequenceName = sequence.name;
		let motifs = sequence.motifs;
		addSequence(genomicContextTable, sequence);		
	}
}

});
