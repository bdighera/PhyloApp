<?php 
include_once('connection.php'); 
$query="select * from student"; 
$result=mysql_query($query); 
?> 
<!DOCTYPE html> 
<html> 
	<head> 
		<title> Fetch Data From Database </title> 
	</head> 
	<body> 
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
		
		<?php while($rows=mysql_fetch_assoc($result)) 
		{ 
		?> 
		<tr> <td><?php echo $rows['ProteinID']; ?></td> 
		<td><?php echo $rows['CommonName']; ?></td> 
		<td><?php echo $rows['ProteinID']; ?></td> 
		<td><?php echo $rows['CDSAccession']; ?></td>
		<td><?php echo $rows['GenomicAccession']; ?></td>
		<td><?php echo $rows['GeneID']; ?></td>
		<td><?php echo $rows['Taxonomy']; ?></td>
		</tr> 
	<?php 
               } 
          ?> 

	</table> 
	</body> 
	</html>