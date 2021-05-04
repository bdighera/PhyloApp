var deletedSequences = [];

var domainHeight = 48;
var domainHeight = 48;
var domainTop = 1;
var domainStart = 26;
var domainEnd = 225;
var domainWidth = domainEnd - domainStart;
var sequenceEnd = 250;
var sequenceStart = 0;
var sequenceTop = 0;
var sequenceBottom = 50;
var sequenceHeight = 20;
var sequenceMidHeight = sequenceHeight/2;
var sequenceWidth = 250;
var intronHeight = 10;
var intronWidth = 10;
var intronWidthOffset = (intronWidth/2)-1;
var intronHeightOffset = (intronHeight/2)-1;

var exportButton = document.getElementById('exportButton');

$.get( "http://localhost:5000/introns", function( data ) {
	$( ".result" ).html( data );
	sequenceData(data);		
});

function postSequences() {
    let postData = '';
    for(let sequence of deletedSequences){
        postData += sequence + ' ';
    }
    alert(deletedSequences);
    $.post("http://localhost:5000/introns", {'deleted_sequences':postData});
}  

function sequenceLine(sequenceContext, end) {
	sequenceContext.strokeStyle = "black";
	sequenceContext.lineWidth = 2;
	sequenceContext.beginPath();
	sequenceContext.moveTo(0, sequenceMidHeight);
	sequenceContext.lineTo(end, sequenceMidHeight);
	sequenceContext.stroke();
}

function addOptions(sequenceRow, sequence) {
	let optionsCell = sequenceRow.insertCell();
	let sequenceNameDiv = document.createElement("div");
	sequenceNameDiv.innerHTML += sequence.name;	
	sequenceNameDiv.className = "sequenceNameStyle";

	let deleteDiv = document.createElement("input");
	deleteDiv.className = "sequenceNameStyle";
	deleteDiv.type = "checkbox";
	deleteDiv.name = sequence.name;
	deleteDiv.addEventListener('change', onClickDelete);

	optionsCell.appendChild(deleteDiv);
	optionsCell.appendChild(sequenceNameDiv);
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

function addIntron(sequenceContext, intron) {
	let intronStart = intron.startLocation;
	let intronEnd = intron.endLocation;
	let intronTop = sequenceMidHeight - 1;

	sequenceContext.fillStyle = intron.background;
	sequenceContext.fillRect(intronStart-intronWidthOffset, intronTop-intronHeightOffset, intronWidth, intronHeight);
}

function addIntrons(sequenceRow, intronSequence) {
	let intronCell = sequenceRow.insertCell();
	let intronDiv = document.createElement("div");
	intronDiv.id = intronSequence.name+"div";
  let intronCanvas = document.createElement('canvas');
	intronCanvas.id = intronSequence.name+"canvas";
  intronCanvas.height = sequenceHeight;
 	intronCanvas.width = sequenceWidth;
	intronDiv.appendChild(intronCanvas);
	intronContext = intronCanvas.getContext('2d');

	let introns = intronSequence.introns;
	let sequenceLength = intronSequence.len;	
  sequenceLine(intronContext, sequenceLength);	
	for (let intron of introns) {
		addIntron(intronContext, intron);
	}
	intronCell.appendChild(intronDiv);
}
function addSequence(intronSequence) {
  var intronTable = document.getElementById('sequenceTable');
  let sequenceRow =  intronTable.insertRow();
	addOptions(sequenceRow, intronSequence);
	addIntrons(sequenceRow, intronSequence);
}


function sequenceData(data) {
	let intronSequences = data.Sequences;
	for (let sequence of intronSequences) {
		addSequence(sequence);		
	}
}
