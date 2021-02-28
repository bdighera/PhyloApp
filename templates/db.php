<?php 
include_once('connection.php'); 
?> 
<!DOCTYPE html> 
<html> 
	<head> 
		<title> Sequence Records </title> 
	</head> 
	<body>
		<?php
		$sql = "SELECT * FROM sequence records;";
		$results = mysqli_query($conn, $sql);
		$resultCheck = mysqli_num_rows($result);
		if ($resultCheck > 0 ) {
			while ($row = mysqli_fetch_assoc()) {
				echo $rows['ProteinID'];
				echo $rows['CommonName'];
				echo $rows['ProteinID'];
				echo $rows['CDSAccession'];
				echo $rows['GenomicAccession'];
				echo $rows['GeneID'];
				echo $rows['Taxonomy'];
			}
		}
		?>
	<table align="center" border="1px" style="width:1000px; line-height:40px;"> 
	<tr> 
		<th colspan="7"><h2> Sequence Records</h2></th> 
		</tr> 
			  <th> Protein Accession </th>
			  <th> Common Name </th>
			  <th> Protein ID </th> 
			  <th> CDS Accession </th>
			  <th> Genomic Accession </th>
			  <th> Gene ID </th>
			  <th> Taxonomy </th>
			  
		</tr> 

	</table> 
	</body> 
	</html>