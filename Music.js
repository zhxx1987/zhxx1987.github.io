
// use global var audio = new Audio();
// otherwise everything will break 
// and you will die, trapped in callback hell, unable to escape, 
// pursued by demons. Lips chapped, you call out for help, 
// but no one can hear you over your whining CPU fan.
// 
// also, make sure to define: 
// function message(msg) {
//   document.getElementById("message").innerHTML = msg;
// }
//
// start everything with this:
// audio.loadSampleAudio();

function Audio(url) {
  this.audioBuffer = null;  // audioContext.createBuffer(request.response, false );
  this.audioContext = null; // window.webkitAudioContext();
  this.source = null;       // audioContext.createBufferSource();
  this.processor = null;    // audioContext.createJavaScriptNode(2048 , 1 , 1 );
  this.analyser = null;     // audioContext.createAnalyser();
  this.playing = false;
  this.fftSize = 128;
  this.audioContext = new window.AudioContext();
}

Audio.prototype.loadSampleAudio = function(url) {
  
  this.source = this.audioContext.createBufferSource();
  this.loadAudioBuffer(url);
}

Audio.prototype.loadAudioBuffer = function (url) {
    var request = new XMLHttpRequest(); // load asynchronously
    request.open("GET", url, true);
    request.responseType = "arraybuffer";
    request.caller = this;
    request.onload = function () {
        var audioData = request.response;
        if (this.caller.audioContext.decodeAudioData) {
            this.caller.audioContext.decodeAudioData(audioData, function (buffer) {
                
                this.caller.source.buffer = buffer;
                this.caller.createAudio();
            }, function (e) {
                console.log(e);
                //message("cannot decode mp3");
            });
        } else {
            var channels = 2;
            var frameCount = this.caller.audioContext.sampleRate * 2.0;

            this.caller.audioBuffer = this.caller.audioContext.createBuffer(2, frameCount, this.caller.audioContext.sampleRate);
            this.caller.source.buffer = this.caller.audioBuffer;
            this.caller.source.loop = true;
            this.caller.source.start(0);
            this.caller.playing = true;
        }
    }
    request.send();
}

Audio.prototype.initAudio = function(data) {
    this.source = this.audioContext.createBufferSource();
   

  if(this.audioContext.decodeAudioData) {
    this.audioContext.decodeAudioData(data, function(buffer){
      this.source.buffer = buffer;
      this.createAudio();
    }, function(e) {
      console.log(e);
      //message("cannot decode mp3");
    });
  } else {
    this.source.buffer = this.audioContext.createBuffer(data, false);
    this.createAudio();
  }
}

Audio.prototype.createAudio = function () {
    //this.processor = this.audioContext.createJavaScriptNode(2048 , 1 , 1 );
    this.analyser = this.audioContext.createAnalyser();
    this.analyser.fftSize = this.fftSize;
    this.source.connect(this.analyser); // Connect audio processing graph
    this.analyser.connect(this.audioContext.destination);
    this.source.connect(this.audioContext.destination);
    this.source.loop = true;
    this.source.noteOn(0);
    this.playing = true;
    //message("OK");
}