<!-- $all = glob('/path/to/dir/*.*'); -->

<!DOCTYPE html>
<html>
<body>

<h1>My first PHP page</h1>

<?php
echo "Hello World!";
$images = glob('https://raw.githubusercontent.com/Shu-Lea-Lai/project-visualmemory-sideexp/master/images/*.jpg');
foreach($images as $image) {echo $image;}

?>

</body>
</html>