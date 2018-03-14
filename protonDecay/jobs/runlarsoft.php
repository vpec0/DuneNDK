<?php
header('Content-Type: text/plain');
header('Content-Disposition: attachment; filename="runlarsoft.sh"');

$nevents=1000;
$run=101;

$mysqli = new mysqli("localhost", "larsoft", "larsoft", "larsoft");
if ($mysqli->correct_errno) {
  echo "Failed to connect to MySQL: (" . $mysqli->connect_errno . ") " . $mysqli->connect_error;
  exit;
}

if (!$mysqli->query("LOCK TABLES larsoft WRITE") ||
    !($res=$mysqli->query("SELECT event from larsoft")))
  {
    echo "Table query failed: (" . $mysqli->errno . ") " . $mysqli->error;
  }
else
  {
    $res->data_seek(0);
    $row=$res->fetch_assoc();
    $data=file_get_contents("runlarsoft.sh");
    $data=str_replace("FIRST", $row['event']*$nevents, $data);
    $data=str_replace("NEVENTS", $nevents, $data);
    $data=str_replace("RUN", $run, $data);
    echo $data;
    if (!$mysqli->query("UPDATE larsoft SET event=".($row['event']+1)))
      { 
	echo "Update failed: (" . $mysqli->errno . ") " . $mysqli->error;
      }
  }
$mysqli->query("UNLOCK TABLES");

?>