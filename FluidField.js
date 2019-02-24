
FluidField = function () {
    this.iterations = 20; 	// solver iterations
    this.dt = 0.10; 		// timestep
    this.dens = null; 	// density field
    this.dens_prev = null;
    this.u = null; 		// velocity x field
    this.u_prev = null;
    this.v = null; 		// velocity y field
    this.v_prev = null;
    this.fx = null;
    this.fy = null;
    this.width = 128; 	// mesh size x
    this.height = 128; 	// mesh size y
    this.rowSize; 		// rowstride
    this.size; 			// mesh points
    this.interpolate = true; // interpolate when setting/getting fields
    this.damp = 1;     // velocity damping
    this.h = 1.0/this.width;
    this.reset();
}

FluidField.prototype.reset = function() {
	this.rowSize = this.width;
	this.size = (this.width) * (this.height);
	this.dens = new Float32Array(this.size);
	this.dens_prev = new Float32Array(this.size);
	this.u = new Float32Array(this.size);
	this.u_prev = new Float32Array(this.size);
	this.v = new Float32Array(this.size);
	this.v_prev = new Float32Array(this.size);
	this.fx = new Float32Array(this.size);
	this.fy = new Float32Array(this.size);
}

FluidField.prototype.setResolution = function (x, y) {
	this.width = x;
	this.height = y;
	this.reset();
}

FluidField.prototype.lerp = function (a, b, alpha) {
    return (1 - alpha) * a + alpha * b;

}

FluidField.prototype.bilerp = function (x, y, field) {

    var xx = x; if (xx >= this.width-1) xx = this.width-1 - 0.01;
    var yy = y; if (yy >= this.height-1) yy = this.height-1 - 0.01;
    if (xx <= 0) xx = 0.01;
    if (yy <= 0) yy = 0.01;
    var x_ = Math.floor(xx);
    var y_ = Math.floor(yy);
    var fx = xx - x_;
    var fy = yy - y_;
    var rowSize = this.rowSize;

    return this.lerp(this.lerp(field[y_ * rowSize + x_], field[y_ * rowSize + x_ + 1], fx),
                this.lerp(field[(y_ + 1) * rowSize + x_], field[(y_ + 1) * rowSize + x_ + 1], fx), fy);

}

FluidField.prototype.bilerpQuad = function (x, y, field, slv, v00, v01, v10, v11) {

    var xx = x; if (xx >= this.width - 1) xx = this.width - 1 - 0.01;
    var yy = y; if (yy >= this.height - 1) yy = this.height - 1 - 0.01;
    if (xx <= 0) xx = 0.01;
    if (yy <= 0) yy = 0.01;
    var x_ = Math.floor(xx);
    var y_ = Math.floor(yy);
    var fx = xx - x_;
    var fy = yy - y_;
    var rowSize = this.rowSize;
    v00 = field[y_ * rowSize + x_];
    v01 = field[y_ * rowSize + x_ + 1];
    v10 = field[(y_ + 1) * rowSize + x_];
    v11 = field[(y_ + 1) * rowSize + x_ + 1];
    slv = this.lerp(this.lerp(v00,v01, fx),
                this.lerp(v10,v11, fx), fy);

}

FluidField.prototype.SL_advection = function (dt, field_new, field_origin, u, v) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    var inv_h = 1.0 / this.h;
    for (var j = 1; j < height-1; j++) {
        for (var i = 1; i < width-1; i++) {
            var px = i + 0.5;
            var py = j + 0.5;
            var vx = u[j * rowSize + i];
            var vy = v[j * rowSize + i];
            var bpx = px - dt * vx * inv_h;
            var bpy = py - dt * vy * inv_h;
            var val = this.bilerp(bpx - 0.5, bpy - 0.5, field_origin);
            field_new[j * rowSize + i] = val;
        }

    }
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            if (i == 0 || j == 0 || i == width - 1 || j == height - 1) {
                var ij = j * rowSize + i;
                field_new[ij] = 0;
            }

        }
    }
}
FluidField.prototype.clampExtrema = function (dt,f_n,f_o,u,v)
{
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    var inv_h = 1.0 / this.h;
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            var px = i + 0.5;
            var py = j + 0.5;
            var vx = u[j * rowSize + i];
            var vy = v[j * rowSize + i];
            var bpx = px - dt * vx * inv_h;
            var bpy = py - dt * vy * inv_h;
            var slv;
            var v00,v01,v10,v11;
            this.bilerpQuad(bpx - 0.5, bpy - 0.5, f_o,slv,v00,v01,v10,v11);
            var min_v = Math.min(Math.min(Math.min(v00,v01),v10),v11);
            var max_v = Math.max(Math.max(Math.max(v00,v01),v10),v11);
            if(f_n[j*rowSize+i]<min_v||f_n[j*rowSize+i]>max_v)
            {
                f_n[j*rowSize+i] = slv;
            }
        }

    }
}
FluidField.prototype.MC_field = function (dt, field_new, field_temp, field_origin, u, v) {
    for (var i = 0; i < this.size; i++) {
        field_temp[i] = 0;
        field_new[i] = 0;
    }
    this.SL_advection(dt, field_temp, field_origin, u, v);
    this.SL_advection(-dt, field_new, field_temp, u, v);
    for (var i = 0; i < this.size; i++) {
        field_temp[i] = field_origin[i] + 0.5 * (field_origin[i] - field_new[i]);
    }
    for (var i = 0; i < this.size; i++) {
        field_new[i] = 0;
    }
    this.SL_advection(dt, field_new, field_temp, u, v);
}
FluidField.prototype.MCAdv = function (dt, u, v, u0, v0) {
    var u_temp = new Float32Array(this.size);
    var v_temp = new Float32Array(this.size);
    this.MC_field(dt, u, u_temp, u0, u0, v0);
    this.MC_field(dt, v, v_temp, v0, u0, v0);
    this.clampExtrema(dt, u, u0, u0, v0);
    this.clampExtrema(dt, v, v0, u0, v0);
}

FluidField.prototype.computeVort = function (w, u, v) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    var inv_h = 1.0 / this.h;
    for (var j = 1; j < height-1; j++) {
        for (var i = 1; i < width-1; i++) {
            var ij = j * rowSize + i;
            w[ij] = (v[ij + 1] - v[ij - 1] - u[ij + rowSize] + u[ij - rowSize]) * 0.5 * inv_h;
        }
    }
}
FluidField.prototype.VortConf = function (w, u, v) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    var inv_h = 1.0 / this.h;
    for (var j = 1; j < height - 1; j++) {
        for (var i = 1; i < width - 1; i++) {
            var ij = j * rowSize + i;
            var gvx = (w[ij + 1] - w[ij - 1]) * 0.5 * inv_h;
            var gvy = (w[ij + rowSize] - w[ij - rowSize]) * 0.5 * inv_h;
            var norm = Math.sqrt(gvx * gvx + gvy * gvy);
            gvx = gvx / (norm + 0.01);
            gvy = gvy / (norm + 0.01);
            var fcx = this.h * 0.15 * (gvy * w[ij]);
            var fcy = this.h * 0.15 * (-gvx * w[ij]);
            u[ij] += fcx;
            v[ij] += fcy;
        }
    }
}
FluidField.prototype.addCurl = function (u, v, psi) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    var inv_h = 1.0 / this.h;
    for (var j = 1; j < height - 1; j++) {
        for (var i = 1; i < width - 1; i++) {
            var ij = j * rowSize + i;
            var ux = 0.5 * inv_h * (psi[ij + rowSize] - psi[ij - rowSize]);
            var uy = 0.5 * inv_h * (psi[ij - 1] - psi[ij + 1]);
            u[ij] += ux;
            v[ij] += uy;
        }
    }
}

FluidField.prototype.addForce = function (u, v, fx, fy, dt) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    for (var i = 0; i < this.size; i++) {
        u[i] += dt * fx[i];
        v[i] += dt * fy[i];
    }
    
}
FluidField.prototype.computeDiv = function (u, v, div) {
    var inv_h = 1.0 / this.h;
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    for (var i = 0; i < this.size; i++) {
        div[i] = 0;
    }

    for (var j = 1; j < height - 1; j++) {
        for (var i = 1; i < width - 1; i++) {
            var left = u[j * rowSize + i - 1];
            var right = u[j * rowSize + i + 1];
            var top = v[(j + 1) * rowSize + i];
            var bot = v[(j - 1) * rowSize + i];

            div[j * rowSize + i] = 0.5 * inv_h * (right - left + top - bot) * this.h * this.h;
        }

    }
}

FluidField.prototype.RBGS = function (x, b, iter, width, height) {
    
    for (var k = 0; k < iter; k++) {
        var rowSize = width;
        for (var j = 0; j < height; j++) {
            for (var i = 0; i < width; i++) {
                if (((i + j) % 2) == 0) {

                var left = (i == 0) ? 0 : x[j * rowSize + i - 1];
                var right = (i == width - 1) ? 0 : x[j * rowSize + i + 1];
                var top = (j == height - 1) ? 0 : x[(j + 1) * rowSize + i];
                var bot = (j == 0) ? 0 : x[(j - 1) * rowSize + i];

                x[j * rowSize + i] = -0.25 * (b[j * rowSize + i] - left - right - top - bot);
                }
            }
        }
        

        for (var j = 0; j < height; j++) {
            for (var i = 0; i < width; i++) {
                if (((i + j) % 2) == 1) {

                    var left = (i == 0) ? 0 : x[j * rowSize + i - 1];
                    var right = (i == width - 1) ? 0 : x[j * rowSize + i + 1];
                    var top = (j == height - 1) ? 0 : x[(j + 1) * rowSize + i];
                    var bot = (j == 0) ? 0 : x[(j - 1) * rowSize + i];

                    x[j * rowSize + i] = -0.25 * (b[j * rowSize + i] - left - right - top - bot);
                }
            }
        }
    }
}
FluidField.prototype.restriction = function (x_c, x_f, width, height)
{
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            var ij = j * width + i;
            var i2 = i * 2;
            var j2 = j * 2;
            var i2j200 = j2 * width * 2 + i2;
            var i2j201 = i2j200 + 1;
            var i2j210 = i2j200 + width * 2;
            var i2j211 = i2j210 + 1;
            x_c[ij] = x_f[i2j200] + x_f[i2j201] + x_f[i2j210] + x_f[i2j211];
        }
    }
}
FluidField.prototype.prolongation = function (x_f, x_c, width, height)
{
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            var ij = j * width + i;
            var i2 = Math.floor(i/2);
            var j2 = Math.floor(j/2);
            x_f[ij] += x_c[i2 + j2 * width / 2];
        }
    }
}
FluidField.prototype.compRes = function (x, b, r, width, height) {


    
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {


            var left = (i == 0) ? 0 : x[j * width + i - 1];
            var right = (i == width - 1) ? 0 : x[j * width + i + 1];
            var top = (j == height - 1) ? 0 : x[(j + 1) * width + i];
            var bot = (j == 0) ? 0 : x[(j - 1) * width + i];

            r[j * width + i] = (b[j * width + i] - left - right - top - bot + 4 * x[j * width + i]);

        }
    }


}
FluidField.prototype.multigrid = function (x, b, width, height) {

    if (width * height <= 64) {
        this.RBGS(x, b, 20, width, height);
    }
    else {
        //RBGS
        this.RBGS(x, b, 2, width, height);
        var r = new Float32Array(width * height);
        this.compRes(x, b, r, width, height);
        var r_c = new Float32Array(width * height / 4);
        this.restriction(r_c, r, width / 2, height / 2);
        var x_c = new Float32Array(width * height / 4);
        this.multigrid(x_c, r_c, width / 2, height / 2);
        this.prolongation(x, x_c, width, height);
        this.RBGS(x, b, 2, width, height);
    }

}


FluidField.prototype.project2 = function (u, v, p, div) {
    var width = this.width;
    var height = this.height;
    var rowSize = this.rowSize;
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            if (i == 0 || j == 0 || i == width - 1 || j == height - 1) {
                var ij = j * rowSize + i;
                u[ij] = 0;
                v[ij] = 0;
            }

        }
    }
    this.computeDiv(u, v, div);
    for (var i = 0; i < this.size; i++) {
        p[i] = 0;
    }
    var inv_h = 1.0 / this.h;
    this.multigrid(p, div, width, height);
    //this.multigrid(p, div, width, height);


    for (var j = 1; j < height - 1; j++) {
        for (var i = 1; i < width - 1; i++) {
            var ij = j * rowSize + i;
            u[ij] = u[ij] - 0.5 * inv_h * (p[ij + 1] - p[ij - 1]);
            v[ij] = v[ij] - 0.5 * inv_h * (p[ij + rowSize] - p[ij - rowSize]);
        }
    }
    for (var j = 0; j < height; j++) {
        for (var i = 0; i < width; i++) {
            if (i == 0 || j == 0 || i == width - 1 || j == height - 1) {
                var ij = j * rowSize + i;
                u[ij] = 0;
                v[ij] = 0;
            }

        }
    }
}




FluidField.prototype.velocityStep = function (u, v, u0, v0, dt, adv) {

    var vort = new Float32Array(this.size);
    
    for (var i = 0; i < this.size; i++) {
        vort[i] = 0;
    }
    this.computeVort(vort, u0, v0);
    if (adv == 0) {

        this.SL_advection(dt, u, u0, u0, v0);
        this.SL_advection(dt, v, v0, u0, v0);
        

    } else {

         this.MCAdv(dt, u, v, u0, v0);
    }
    this.addForce(u, v, this.fx, this.fy, dt);
    this.VortConf(vort, u, v);
    this.project2(u, v, u0, v0);

    for (var i = 0; i < this.size; i++) {
        u0[i] = u[i];
        v0[i] = v[i];
    }
    for (var i = 0; i < this.size; i++) {
        u0[i] *= 0.985;
        v0[i] *= 0.985;
        u[i] *= 0.985;
        v[i] *= 0.985;
    }
}


FluidField.prototype.setVelocity = function (x, y, xv, yv) {

    for (var j = y - 1; j <= y + 1; j++) {
        for (var i = x - 1; i <= x + 1; i++) {
            if (i > 0 && i < this.width - 1 && j > 0 && j < this.height - 1) {
                this.fx[(i) + (j) * this.rowSize] += xv;
                this.fy[(i) + (j) * this.rowSize] += yv;
            }
        }

    }



}

FluidField.prototype.addBouyancy = function (x, y, xv, yv) {

    var i = Math.floor(x - 0.5);
    var j = Math.floor(y - 0.5); 
    if(i>0&&i<this.width&&j>0&&j<this.height)
    {

        this.fx[i + j * this.rowSize] += xv;
        this.fy[i + j * this.rowSize] += yv;
    }




}

FluidField.prototype.getXVelocity = function (x, y) {
    return this.bilerp(x, y, this.u);
}

FluidField.prototype.getYVelocity = function (x, y) {
    return this.bilerp(x, y, this.v);
}

FluidField.prototype.update = function () {
    if (this.dt <= 0.05) {

        this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, this.dt, 1);
        for (var i = 0; i < this.size; i++) {
            this.fx[i] = 0;
            this.fy[i] = 0;
        }
    }
    else if (this.dt <= 0.1) {
        this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, this.dt,0);
        for (var i = 0; i < this.size; i++) {
            this.fx[i] = 0;
            this.fy[i] = 0;
        }
    }
    else {
        this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, 0.1, 1);
        this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, this.dt - 0.1, 0);
        for (var i = 0; i < this.size; i++) {
            this.fx[i] = 0;
            this.fy[i] = 0;
        }
    }
    //this.velocityStep(this.u, this.v, this.u_prev, this.v_prev, 0.5 * this.dt);
}