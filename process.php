<?php
session_start();
if( isset($_POST) ){
     
    //form validation vars
    $formok = true;
    $errors = array(); 
   
    
     
    //form data
    $name = $_POST['name'];    
    $email = $_POST['email'];
    $telephone = $_POST['telephone'];
    $arrivaltime = $_POST['arrivaltime'];
    $arrivaldate = $_POST['arrivaldate'];
	  $invitecode = $_POST['invitecode'];
    $wechatID   = $_POST['wechatID'];
    $flightnum  = $_POST['flightnum'];
    
	if(empty($name)){
        $formok = false;
        echo "you have not entered a name";
    }
     
    //validate email address is not empty
    if(empty($email)){
        $formok = false;
        echo 'you have not entered an email';
    //validate email address is valid
    }elseif(!filter_var($email, FILTER_VALIDATE_EMAIL)){
        $formok = false;
        echo "the email address is not valid";
    }
    
    
    list($y, $m, $d) = explode('-', $arrivaldate);
    
    if(empty($arrivaldate))
    {
        $formok = false;
        echo "you have not entered a date.<bg>";
    }elseif($y!=2017)
    {
        $formok = false;
        echo "wrong year!.<bg>";
    }elseif($m!=9)
    {
        $formok = false;
        echo "wrong month!.<bg>";
    }elseif($d>17)
    {
        $formok = false;
        echo "you are going to miss the event! .<bg>";
    }
    
    


	if($formok){
        echo "success! an email copy has been sent to both you and me! see you in Bali!";
         
        $emailbody = "You have recieved a new schedule, 
                      you are going to pick up  {$name} at {$arrivaldate} on {$arrivaltime} UTC+08:00,
                      his/her Email Address is {$email}
                      his/her Telephone is  {$telephone}";
					  
        $emailbody2 = "Xinxin have recieved your schedule, 
                      you are going to meet him at {$arrivaldate} on {$arrivaltime} UTC+08:00,
                      you can reach him by 7788929102"; 
                      
        mail("zhangshinshin@gmail.com","pickup schedule", $emailbody);
        mail($email,"pickup schedule", $emailbody2);
         
    }
 
  if($invitecode == "0xaed55")
  {
    $myfile = fopen("schedule.txt", "a") or die("Unable to open file!");
    $txt = sprintf("%-30s", $name);
    fwrite($myfile, $txt);
    $txt = sprintf("%-30s", "{$arrivaldate} {$arrivaltime} UTC+08:00");
    fwrite($myfile, $txt);
    $txt = sprintf("%-15s", $flightnum);
    fwrite($myfile, $txt);
    $txt = sprintf("%-20s", $telephone);
    fwrite($myfile, $txt);
    $txt = sprintf("%-20s", $wechatID);
    fwrite($myfile, $txt);
    fwrite($myfile, "\r\n");
    fclose($myfile);
   }
	  
}
?>