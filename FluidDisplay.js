
function FluidDisplay(field) {
  this.additive = true;
  this.canvas = document.getElementById("canvas-fluid");
  this.context = this.canvas.getContext("2d");
  this.image = null; // created in renderParticles
  // buffer is for displaying density
  this.buffer = document.createElement("canvas");
  this.buffer.width = field.width;
  this.buffer.height = field.height;
  this.buffercontext = this.buffer.getContext("2d");
  this.bufferimage = this.buffercontext.createImageData(field.width, field.height);
  // set alpha to 255
  
    for (var i = 3; i < field.width * field.height * 4; i += 4) {
        this.bufferimage.data[i] = 255;
    }
}

FluidDisplay.prototype.clear = function() {
  this.context.fillStyle = 'black';
  this.context.fillRect(0, 0, this.canvas.width, this.canvas.height);
}

FluidDisplay.prototype.renderDensityField = function(field) {
  for (var x = 0; x < field.width; x++) {
    for (var y = 0; y < field.height; y++) {
      var d = field.getDensity(x, y) / 5;
      if (d > 1) d = 1;
      var base = 4 * (y * field.width + x);
      d = d * 2 - 1;
      this.bufferimage.data[base + 0] = jet.r(d);
      this.bufferimage.data[base + 1] = jet.g(d);
      this.bufferimage.data[base + 2] = jet.b(d);
    }
  }
  this.buffercontext.putImageData(this.bufferimage, 0, 0);
  var hscale = this.canvas.width / this.buffer.width;
  var vscale = this.canvas.height / this.buffer.height;
  this.context.scale(hscale, vscale);
  this.context.drawImage(this.buffer, 0, 0);
  this.context.scale(1 / hscale, 1 / vscale);
}

FluidDisplay.prototype.renderVelocityField = function(field) {
  this.context.save();
  this.context.lineWidth = 1;
  var wscale = this.canvas.width / field.width;
  var hscale = this.canvas.height / field.height;
  this.context.strokeStyle = 'white';
  var scale = 20;
  this.context.beginPath();
  for (var x = 0; x < field.width; x++) {
    for (var y = 0; y < field.height; y++) {
      var vx = field.getXVelocity(x, y);
      var vy = field.getYVelocity(x, y);
      this.context.moveTo(x * wscale + 0.5 * wscale, y * hscale + 0.5 * hscale);
      this.context.lineTo((x + 0.5 + scale * vx) * wscale,
                          (y + 0.5 + scale * vy) * hscale);
    }
  }
  this.context.stroke();
  this.context.restore();
}

FluidDisplay.prototype.renderParticles = function (field, light, px, py, pc, pl) {

    var coef_lookup = new Float32Array(9);
    for (var jj = -1; jj <= 1; jj++) {
        for (var ii = -1; ii <= 1; ii++) {
            var j = jj + 1;
            var i = ii + 1;
            coef_lookup[i + j * 3] = Math.exp(-(ii * ii - jj * jj) * 0.00001);
        }
    }
    this.image = this.context.createImageData(this.canvas.width, this.canvas.height);
    var wscale = this.canvas.width / field.width;
    var hscale = this.canvas.height / field.height;
    var n = px.length;
    var h, s, v, r, g, b, j, f, p, q, t;
    for (var i = 0; i < n; i++) {
        var cxx = Math.floor(px[i] * wscale);
        var cyy = Math.floor(py[i] * hscale);
        //for (var jj = -1; jj <= 1; jj++)
        //for (var ii = -1; ii <= 1; ii++) {
        var cx = cxx; // ;
        var cy = cyy; // +jj;
        cx = Math.min(Math.max(0, cx), this.canvas.width);
        cy = Math.min(Math.max(0, cy), this.canvas.width);
        var coef = 1;// coef_lookup[ii + 1 + (jj + 1) * 3];

        var base = 4 * (cy * this.canvas.width + cx);
        h = pc[i];
        s = pl[i] / life;
        v = light * 2 * coef * s * s;
        { // hsv -> rgb, inlined for performance
            j = Math.floor(h * 6);
            f = h * 6 - j;
            p = v * (1 - s);
            q = v * (1 - f * s);
            t = v * (1 - (1 - f) * s);
            switch (j % 6) {
                case 0: r = v, g = t, b = p; break;
                case 1: r = q, g = v, b = p; break;
                case 2: r = p, g = v, b = t; break;
                case 3: r = p, g = q, b = v; break;
                case 4: r = t, g = p, b = v; break;
                case 5: r = v, g = p, b = q; break;
            }
        }
        // alpha
        this.image.data[base + 3] = 255;
        // blending modes
        if (this.additive) {
            this.image.data[base] += Math.floor(r * 255);
            this.image.data[base + 1] += Math.floor(g * 255);
            this.image.data[base + 2] += Math.floor(b * 255);
        } else {
            this.image.data[base] = Math.floor(r * 255);
            this.image.data[base + 1] = Math.floor(g * 255);
            this.image.data[base + 2] = Math.floor(b * 255);
        }
        //}
    }
    this.context.putImageData(this.image, 0, 0);
}