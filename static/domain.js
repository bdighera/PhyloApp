var deletedSequences = [];

$.get( "http://localhost:5000/domains", function( data ) {
    $( ".result" ).html( data );

    sequenceData(data);

});

function postSequences() {
    let postData = '';
    for(let sequence of deletedSequences){
        postData += sequence + ' ';
    }
    alert(deletedSequences);
    $.post("http://localhost:5000/domains", {'deleted_sequences':postData});
}

function sequenceData(data) {
    var names = "";
    let sequences = data.Sequences;
    var seqTable = document.getElementById('sequenceTable');
    var maxWidth = 0;
    for(let sequence of sequences) {
      for(let domain of sequence.domains) {
        if(domain.endLocation >= maxWidth) { 
          maxWidth = domain.endLocation;
        }
      }
    }
    //var deletedSequences = [];

    for(let sequence of sequences){
        names = names + sequence.name + " ";

        let seqRow = seqTable.insertRow();
        let deleteCell = seqRow.insertCell();
        let nameCell = seqRow.insertCell();
        let domainCell = seqRow.insertCell();

        //var midHeight = sequence.domains[0].height/2;
        //var midWidth = sequence.domains[0].startLocation + ((sequence.domains[0].endLocation - sequence.domains[0].startLocation)/2);
        //var startLocation = sequence.domains[0].startLocation;
        //var endLocation = sequence.domains[0].endLocation;
        //var height = sequence.domains[0].height;

        let deleteCheck = document.createElement('input');
        deleteCheck.type = "checkbox";
        deleteCheck.value = sequence.name;
        deleteCheck.id = sequence.name + 'check';
        deleteCheck.addEventListener('change', function() {
            if(this.checked) {
                deletedSequences.push(deleteCheck.value);
                alert(deletedSequences);
                document.getElementById(deleteCheck.value + 'canvas').style.opacity = "0.25";
            } else {
                let deletedLocation = deletedSequences.indexOf(deleteCheck.value);
                if (deletedLocation > -1) {
                    deletedSequences.splice(deletedLocation, 1);
                }
                alert(deletedSequences);
                document.getElementById(deleteCheck.value + 'canvas').style.opacity = "1.00";
            }
        });

        deleteCell.appendChild(deleteCheck);


        let nameDiv = document.createElement('div');
        nameDiv.id = sequence.name;
        nameDiv.innerHTML = sequence.name;
        nameDiv.className = "sequenceNameStyle";
        nameCell.appendChild(nameDiv);


        let nameCanvas = document.createElement('canvas');
        nameCanvas.id = sequence.name + 'canvas';
        nameCanvas.width = maxWidth;
        nameCanvas.height = sequence.domains[0].height;
        nameCanvas.className = "sequenceCanvasStyle";
        domainCell.appendChild(nameCanvas);

        //let canvas = document.getElementById(nameCanvas.id);
        if (nameCanvas.getContext) {
            var ctx = nameCanvas.getContext('2d');
            for(let domain of sequence.domains) {
              let midHeight = domain.height/2;
              let startLocation = domain.startLocation;
              let endLocation = domain.endLocation;
              let midWidth = startLocation + ((endLocation - startLocation)/2);
              let height = domain.height;

              ctx.beginPath();
              ctx.moveTo(startLocation, midHeight);
              ctx.lineTo(midWidth, 0);
              ctx.lineTo(endLocation, midHeight);
              ctx.lineTo(midWidth, height);
              ctx.fillStyle=domain.background;
              ctx.fill();
           }
        }
    } 
}
