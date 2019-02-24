
function BeatDetector() {
  this.history = [];
  this.historyBins = 12;
  this.beatChance = 0;
  this.binChance = [];
  this.historyWindow = 30;
  this.levelHistory = [];
  this.threshold = 2; // times stdev 
  for (var i = 0; i < this.historyBins; i++) {
    this.history[i] = []; 
    this.beatChance[i] = 0;
  }
}

BeatDetector.prototype.sample = function(analyser) {
  var fftSize = analyser.frequencyBinCount;
  var freqByteData = new Uint8Array(fftSize);
  analyser.smoothingTimeConstant = 0.1;
  analyser.getByteFrequencyData(freqByteData);
  // get average level
  var sum = 0;
  for (var i = 0; i < fftSize; i++) {
    var weight = 2 - (i / fftSize);
    sum += weight * freqByteData[i];
  }
  var average = sum / fftSize / 255;
  this.levelHistory.push(average);  
  if (this.levelHistory.length > this.historyWindow) {
    this.levelHistory.shift();
  }
  // add bins to history
  var binsPerBin = Math.floor(fftSize / this.historyBins);
  for (var i = 0; i < this.historyBins; i++) {
    var sum = 0;
    for (var j = 0; j < binsPerBin; j++) {
      sum += freqByteData[i * binsPerBin + j] / 255;
    }
    this.history[i].push(sum / binsPerBin);
    if (this.history[i].length > this.historyWindow) {
      this.history[i].shift();
    }
  }
  // check if newest entry in history is significant
  if (this.levelHistory.length >= this.historyWindow) {
    
    // first look at general levels
    var cur = this.levelHistory[this.historyWindow - 1];
    var mean = 0;
    for (var i = 0; i < this.historyWindow - 1; i++) {
      mean += this.levelHistory[i];
    }
    mean /= this.levelHistory.length; 
    var variance = 0;
    for (var i = 0; i < this.historyWindow - 1; i++) {
      variance += (this.levelHistory[i] - mean) * (this.levelHistory[i] - mean);
    }
    var stdev = Math.sqrt(variance / this.historyWindow);
    this.beatChance = Math.abs(cur - mean) / (stdev + 0.01);
    
    // then look at each bin
    for (var i = 0; i < this.historyBins; i++) {
      var cur = this.history[i][this.historyWindow - 1];
      var mean = 0;
      for (var j = 0; j < this.historyWindow - 1; j++) {
        mean += this.history[i][j];
      }
      mean /= this.history[i].length;
      var variance = 0; // sample variance
      for (var j = 0; j < this.historyWindow - 1; j++) {
        variance += (this.history[i][j] - mean) * (this.history[i][j] - mean);
      }
      variance /= this.historyWindow;
      var stdev = Math.sqrt(variance);
      var weight = 2 - (i / this.historyBins);
      this.binChance[i] = Math.abs(cur - mean) / (stdev + 0.01);
    }
    //console.log(this.beatChance);
  }
}