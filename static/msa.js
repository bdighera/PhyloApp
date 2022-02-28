/**
 * Created by Bryan on 2/27/22.
 */

function formatMSA(msa) {

    seqArray = Array.from(msa)
    // console.log(seqArray)
    seqArray.forEach(function (item, index) {
        const square = document.createElement("div");
        const newContent = document.createTextNode(item);
        square.appendChild(newContent);
        // console.log(item, index);
        const currentDiv = document.getElementById("msa_subset");
        // document.body.appendChild(square);
    });
}

