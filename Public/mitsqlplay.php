<?

$dbuser="r_chen";
$dbpass="mitsqlpass";
$dbname="r_chen+STEM"; #the name of the database
$chandle = mysql_connect("sql.mit.edu", $dbuser, $dbpass)
    or die("Connection Failure to Database");
echo "Connected to database server<br>";
mysql_select_db($dbname, $chandle) or die ($dbname . " Database not found." . $dbuser);
echo "Database " .  $dbname . " is selected<br>";

$var = "sample ";
$trimmed = trim($var); //trim whitespace from the stored variable

// rows to return
$limit=10; 

// check for an empty string and display a message.
if ($trimmed == "")
  {
  echo "<p>Please enter a search...</p>";
  exit;
  }

// check for a search parameter
if (!isset($var))
  {
  echo "<p>We dont seem to have a search parameter!</p>";
  exit;
  }

$table = "videos";
$field = "keyword";
$resultField = "url";
$resultField1 = 
###########################EDIT SQL QUERY
// Build SQL Query  
$query = "select * from $table where $field like \"%$trimmed%\"  
  order by $field"; // EDIT HERE and specify your table and field names for the SQL query

 $numresults=mysql_query($query);
 $numrows=mysql_num_rows($numresults);

##echo "FOUND ".$numresults." results!";


// next determine if s has been passed to script, if not use 0
  if (empty($s)) {
  $s=0;
  }

// get results
  ####$query .= " limit $s,$limit";##dont set a limit on number of rows to return
  $result = mysql_query($query) or die("Couldn't execute query");

// display what the person searched for
echo "<p>You searched for: &quot;" . $var . "&quot;</p>";

// begin to show results set
echo "Results<br><br>";
$count = 1 + $s ;

// now you can display the results returned//here, display a link to the video page, and add a description under the link, with value in $row["shortdescrip"]; //ie., the "short description" in our VIDEOS table
  while ($row= mysql_fetch_array($result)) {
  $title = $row["$resultField"];
  $title1 = $row["embedtag"];
  $title2 = $row["thumbnail"];
  $title3 = $row["shortdescrip"];
  $title4 = $row["randomfield"];

  echo "$count.)&nbsp;$title<br>$title1<br>$title2<br>$title3<br>$title4<br><br>" ;
  $count++ ;
  }

echo "-------------------------------------------------<br>";

###########OMIT THE FOLLOWING IF WE WANT ALL SEARCH RESULTS TO APPEAR ON ONE PAGE
$currPage = (($s/$limit) + 1);

//break before paging
  echo "<br />";

  // next we need to do the links to other results
  if ($s>=1) { // bypass PREV link if s is 0
  $prevs=($s-$limit);
  print "&nbsp;<a href=\"$PHP_SELF?s=$prevs&q=$var\">&lt;&lt; 
  Prev 10</a>&nbsp&nbsp;";
  }

// calculate number of pages needing links
  $pages=intval($numrows/$limit);

// $pages now contains int of pages needed unless there is a remainder from division

  if ($numrows%$limit) {
  // has remainder so add one page
  $pages++;
  }

// check to see if last page
  if (!((($s+$limit)/$limit)==$pages) && $pages!=1) {

  // not last page so give NEXT link
  $news=$s+$limit;

  echo "&nbsp;<a href=\"$PHP_SELF?s=$news&q=$var\">Next 10 &gt;&gt;</a>";
  }

$a = $s + ($limit) ;
  if ($a > $numrows) { $a = $numrows ; }
  $b = $s + 1 ;
  echo "<p>Showing results $b to $a of $numrows</p>";


mysql_close($chandle);


?>
