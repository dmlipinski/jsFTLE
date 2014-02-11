/*    
@licstart  The following is the entire license notice for the 
JavaScript code in this page.


Copyright (c) 2014, Doug Lipinski, dmlipinski@gmail.com
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


@licend  The above is the entire license notice
for the JavaScript code in this page.
*/


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
|                                                              |
|                       object definitions                     |
|                                                              |
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

function Points(param)
{
	this.param = param;
	this.x = new Float32Array(this.param.N);
	this.y = new Float32Array(this.param.N);

	/*~~~~~~~~~~*\
	|  METHODS:  |
	\*~~~~~~~~~~*/
	this.Clone = Clone;
	function Clone(src)
	{
		for (var i=0; i<this.param.N; i++)
		{
			this.x[i]=src.x[i];
			this.y[i]=src.y[i];
		}
	}

	this.initialize = initialize;
	function initialize()
	{
		var index=0;
		for (var j=0; j<this.param.ny; j++) {
			for (var i=0; i<this.param.nx; i++) {
				this.x[index]=this.param.xmin+i*this.param.dx;
				this.y[index++]=this.param.ymin+j*this.param.dy;
			}
		}
	}

	this.draw = draw;
	function draw(canvas,context,drawingParam,percent) {
		//reset canvas
		context.clearRect(0,0,canvas.width,canvas.height);
		
		context.fillStyle = "#000";
		for (var i=0; i<this.param.N; i++) {
			context.fillRect(
				Math.round(drawingParam.stretch*(this.x[i]-this.param.xmin)-drawingParam.ptsize/2),
				Math.round(canvas.height-drawingParam.stretch*(this.y[i]-this.param.ymin)-drawingParam.ptsize/2),
				drawingParam.ptsize, drawingParam.ptsize);
		}

		//progress bar
		context.fillStyle = "rgba(255,0,0,.5)";
		context.fillRect(0,canvas.height-5,percent*canvas.width,5)
	}
}

function velField(points)
{
	this.param = new Object();
	this.param.nx = points.param.nx;
	this.param.ny = points.param.ny;
	this.param.N = points.param.N;
	this.u = new Float32Array(this.param.N);
	this.v = new Float32Array(this.param.N);

	this.update = update;
	function update(points,t)
	{
		var uv=[0,0]
		for (i=0; i<this.param.N; i++)
		{
			uv = velocity(points.x[i],points.y[i],t);
			this.u[i]=uv[0];
			this.v[i]=uv[1];
		}
	}
}

function ftleField(points)
{
	this.param = new Object();
	this.param.nx = points.param.nx;
	this.param.ny = points.param.ny;
	this.param.N = points.param.N;
	this.val = new Float32Array(this.param.N);
	this.maxFTLE = 0;

	this.reset = reset;
	function reset() {
		this.maxFTLE = 0;
	}

	//compute FTLE values
	this.computeFTLE = computeFTLE;
	function computeFTLE(state) {

		var d11, d12, d21, d22, A, B, C;
		var twoDx=2*state.pointParam.dx, twoDy=2*state.pointParam.dy;

		this.reset();
		var index;
		for (var j=1; j<state.pointParam.ny-1; j++)
		{
			for (var i=1; i<state.pointParam.nx-1; i++)
			{
				index = i+state.pointParam.nx*j;
				d11 = ( state.points.x[index+1]-state.points.x[index-1] )/twoDx;
				d12 = ( state.points.x[index+state.pointParam.nx]-state.points.x[index-state.pointParam.nx] )/twoDy;
				d21 = ( state.points.y[index+1]-state.points.y[index-1] )/twoDx;
				d22 = ( state.points.y[index+state.pointParam.nx]-state.points.y[index-state.pointParam.nx] )/twoDy;

				A = d11*d11+d12*d12;
				B = d11*d21+d12*d22;
				C = d21*d21+d22*d22;

				//do not take the natural log or divide by 2*T, this will be accounted for in the colormap
				this.val[index]=(A+C+Math.sqrt(Math.pow(A+C,2)-4*(A*C-B*B)))/2;
				if (this.val[index]>this.maxFTLE) {this.maxFTLE=this.val[index];}
			}
		}
	}
}

function State(pointParam,timeParam)
{
	this.pointParam=pointParam;
	this.timeParam=timeParam;
	this.points=new Points(pointParam);
	this.vel=new velField(this.points);
	this.fwdftle=new ftleField(this.points);
	this.bkwdftle=new ftleField(this.points);

	this.t=timeParam.t0;
	this.tIndex=0; //current time index
	this.dir=1; //forward (1) or backward (-1) advection
	this.nSteps=Math.round(timeParam.T/timeParam.dt); //number of time steps to take for advection
	this.dtActual=timeParam.T/this.nSteps;

	this.setup=setup;
	function setup()
	{
		this.dir=1;
		this.t=this.timeParam.t0;
		this.tIndex=0;
	}

	//step the simulation forward (or backward, depending on "this.dir") by "dt"
	this.Step = Step;
	function Step()
	{
		var t = this.t;
		var dt=this.dtActual*this.dir;

		//RK4 coefficients:
		var k1=new velField(this.points);
		var k2=new velField(this.points);
		var k3=new velField(this.points);
		var k4=new velField(this.points);

		//var pts = this.points;
		var pts = new Points(this.points.param);
		pts.Clone(this.points);
		k1.update(pts,t);

		eulerStep(pts,k1,dt/2);
		k2.update(pts,t+dt/2)

		pts.Clone(this.points);
		eulerStep(pts,k2,dt/2);
		k3.update(pts,t+dt/2)

		pts.Clone(this.points);
		eulerStep(pts,k3,dt);
		k4.update(pts,t+dt)

		//advect the points:
		RK4Step(this.points,k1,k2,k3,k4,dt);

		//update t
		this.t+=dt;
		this.tIndex++;
	}

	function eulerStep(pts,vel,dt)
	{
		//euler step:
		for (var i=0; i<pts.param.N; i++)
		{
			pts.x[i]=pts.x[i]+dt*vel.u[i];
			pts.y[i]=pts.y[i]+dt*vel.v[i];
		}

	}

	function RK4Step(pts,k1,k2,k3,k4,dt)
	{
		//RK4 step:
		var delta = dt/6;
		for (var i=0; i<pts.param.N; i++)
		{
			//pts.x[i]=pts.x[i]+delta*k1.u[i];
			//pts.y[i]=pts.y[i]+delta*k1.v[i];
			pts.x[i]=pts.x[i]+delta*(k1.u[i]+2*k2.u[i]+2*k3.u[i]+k4.u[i]);
			pts.y[i]=pts.y[i]+delta*(k1.v[i]+2*k2.v[i]+2*k3.v[i]+k4.v[i]);
		}

	}
}





/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
|                                                              |
|                         start program                        |
|                                                              |
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

var myInterval1, myInterval2, running=false;

var velParam = {
	A:0.1,
	eps:0.25,
	omega:Math.PI*2/10
};

//ftle parameters:
var pointParam = {
	nx:201,
	ny:101,
	xmin:0,
	xmax:2,
	ymin:0,
	ymax:1
}
var timeParam = {
	t0:0,
	T:12,
	dt:1
};

document.getElementById("nx").value = pointParam.nx.toString();
document.getElementById("ny").value = pointParam.ny.toString();
document.getElementById("t0").value = timeParam.t0.toString();
document.getElementById("T").value = timeParam.T.toString();
document.getElementById("dt").value = timeParam.dt.toString();






//terms for coordinate change when drawing to screen coordinates:
var drawingParam = {
	ptsize:2, //point size in pixels
	maxWidth:628, //canvas max width
	maxHeight:480 //canvas max height
}


//get drawing context:
var canvas = document.getElementById("main_canvas");
var ctx=canvas.getContext('2d');

//scale canvas to correct aspect ratio
if ( (pointParam.xmax-pointParam.xmin)/(pointParam.ymax-pointParam.ymin)>=drawingParam.maxWidth/drawingParam.maxHeight )
{
	canvas.width=drawingParam.maxWidth;
	canvas.height=drawingParam.maxWidth*(pointParam.ymax-pointParam.ymin)/(pointParam.xmax-pointParam.xmin);
	drawingParam.stretch=drawingParam.maxWidth/(pointParam.xmax-pointParam.xmin);
} else
{
	canvas.height=drawingparam.maxHeight;
	canvas.width=drawingParam.maxHeight*(pointParam.xmax-pointParam.xmin)/(pointParam.ymax-pointParam.ymin);
	drawingParam.stretch=drawingParam.maxHeight/(pointParam.ymax-pointParam.ymin);
}


//MAIN FUNCTION, COMPUTE LCS

function computeLCS() {

	//if the function is already running, kill it before starting again.
	if (running)
	{
		clearInterval(myInterval1);
		clearInterval(myInterval2);
	}
	running=true;

	var test = getParameters();
	if (test!=0) return;

	//derived parameters:
	pointParam.dx=(pointParam.xmax-pointParam.xmin)/(pointParam.nx-1);
	pointParam.dy=(pointParam.ymax-pointParam.ymin)/(pointParam.ny-1);
	pointParam.N=pointParam.nx*pointParam.ny;

	document.getElementById("status").innerHTML = "Performing particle advections...forward time";
	document.getElementById("status").style.color = "black";

	//create objects 
	var state = new State(pointParam,timeParam);

	//forward FTLE computations:
	state.points.initialize();
	myInterval1 = setInterval(
		function() {
			state.Step();
			state.points.draw(canvas,ctx,drawingParam,(state.t-state.timeParam.t0)/state.timeParam.T/2);
			if (state.tIndex==state.nSteps) 
			{
				state.fwdftle.computeFTLE(state);
				clearInterval(myInterval1);

	//backward FTLE computations:
	document.getElementById("status").innerHTML = "Performing particle advections...backward time";
	document.getElementById("status").style.color = "black";
	state.setup();
	state.dir=-1;
	state.points.initialize();
	myInterval2 = setInterval(
		function() {
			state.Step();
			state.points.draw(canvas,ctx,drawingParam,(state.timeParam.t0-state.t)/state.timeParam.T/2+.5);
			if (state.tIndex==state.nSteps) 
			{
				state.bkwdftle.computeFTLE(state);
				clearInterval(myInterval2);
				drawFTLE(state);
				document.getElementById("status").innerHTML = "Done, ready.";
				document.getElementById("status").style.color = "black";
				running=false;
			}
		}
		,1000/60);
			}
		}
		,1000/60);

}



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*\
|                                                              |
|                 start function definitions                   |
|                                                              |
\*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//velocity function
function velocity(x,y,t) {

	var a = velParam.eps*Math.sin(velParam.omega*t),
	b = 1-2*velParam.eps*Math.sin(velParam.omega*t);

	var f, dfdx;

	f = a*x*x+b*x;
	dfdx = 2*a*x+b;

	var uv = [0, 0];
	uv[0] = -Math.PI*velParam.A*Math.sin(Math.PI*f)*Math.cos(Math.PI*y),
	uv[1] = Math.PI*velParam.A*Math.cos(Math.PI*f)*Math.sin(Math.PI*y)*dfdx;

	return uv;
}

















function drawFTLE(state) {

	//to avoid aliasing issues, create the image at 1:1 pixel resolution then scale it at the end
	var W=canvas.width;
	var H=canvas.height;
	canvas.width=state.pointParam.nx;
	canvas.height=state.pointParam.ny;

	//clear canvas
	ctx.clearRect(0,0,state.pointParam.nx+1,state.pointParam.ny+1);

	var index, red, blue;

	//plot FTLE values
	for (var j=1; j<state.pointParam.ny-1; j++)
	{
		for (var i=1; i<state.pointParam.nx-1; i++)
		{
			index = i+state.pointParam.nx*j;
			
			//color determined by the ftle values
			red = Math.log(state.bkwdftle.val[index])/Math.log(state.bkwdftle.maxFTLE);
			blue = Math.log(state.fwdftle.val[index])/Math.log(state.fwdftle.maxFTLE);

			//increase contrast between low and high values by squaring the rg values
			red*=red;
			blue*=blue;

			//use bilinear interpolation for a 2d colormap based on red and blue components
			//red=blue=0 => white, red=blue=1 => purple=rgb(.7,0,.9)
			//rgb = (1-blue)*(red*[1, 0, 0]+(1-red)*[1, 1, 1])+blue*(red*[.7, 0, .9]+(1-red)*[0, 0, 1]))
			ctx.fillStyle = "rgb(".concat(
					Math.floor(256*( (1-blue)*red*1 + (1-blue)*(1-red)*1 + blue*red*.6 )).toString(),",", 
					Math.floor(256*( (1-blue)*(1-red)*1 )).toString(),",",
					Math.floor(256*( (1-blue)*(1-red)*1 + blue*red*.9 + blue*(1-red)*1 )).toString(),")");

			//draw a 1 pixel rectangle
			ctx.fillRect(i,state.pointParam.ny-1-j,1,1);
		}
	}

	//image the canvas
	var image = new Image();
	image.src = canvas.toDataURL("image/png");

	//rescale the canvas and draw the image
	canvas.width=W;
	canvas.height=H;
	ctx.drawImage(image,0,0,W,H);

	//Create a new <img> in the results <div>
	var resultsDiv=document.getElementById('old_results');
	resultsDiv.insertAdjacentHTML('afterbegin','<hr><p>Parameters: A='+velParam.A+', &epsilon;='+velParam.eps+', &omega;='+velParam.omega+', nx='+pointParam.nx+', ny='+pointParam.ny+', t<sub>0</sub>='+timeParam.t0+', T='+timeParam.T+', &Delta;t='+timeParam.dt+'<br><img class="results"></p>');
	var newImg=document.getElementsByClassName("results");
	newImg[0].src=image.src;
	newImg[0].width=W*.75;
	newImg[0].height=H*.75;

}

function getParameters()
{
	//Get inputs from the text boxes and validate before returning.
	// return value of 0 indicates success. return value of 1 indicates validation failure.

	velParam.A = parseFloat(document.getElementById("A").value);
	if (isNaN(velParam.A))
	{
		document.getElementById("status").innerHTML = "ERROR, 'A' must be a number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	velParam.eps = parseFloat(document.getElementById("epsilon").value);
	if (isNaN(velParam.eps) || velParam.eps<0)
	{
		document.getElementById("status").innerHTML = "ERROR, '&epsilon;' must be a non-negative number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	velParam.omega = parseFloat(document.getElementById("omega").value);
	if (isNaN(velParam.omega))
	{
		document.getElementById("status").innerHTML = "ERROR, '&omega;' must be a number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	pointParam.nx = parseInt(document.getElementById("nx").value);
	if (isNaN(pointParam.nx) || pointParam.nx<3)
	{
		document.getElementById("status").innerHTML = "ERROR, 'nx' must be an integer greater than 2.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	pointParam.ny = parseInt(document.getElementById("ny").value);
	if (isNaN(pointParam.ny) || pointParam.ny<3)
	{
		document.getElementById("status").innerHTML = "ERROR, 'ny' must be an integer greater than 2.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	timeParam.t0 = parseFloat(document.getElementById("t0").value);
	if (isNaN(timeParam.t0))
	{
		document.getElementById("status").innerHTML = "ERROR, 't<sub>0</sub>' must be a number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}

	timeParam.T = parseFloat(document.getElementById("T").value);
	if (isNaN(timeParam.T) || timeParam.T<=0)
	{
		document.getElementById("status").innerHTML = "ERROR, 'T' must be a positive number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}
	
	timeParam.dt = parseFloat(document.getElementById("dt").value);
	if (isNaN(timeParam.dt) || timeParam.dt<=0)
	{
		document.getElementById("status").innerHTML = "ERROR, '&Delta;t' must be a positive number.";
		document.getElementById("status").style.color = "red";
		return 1;
	}


	return 0;
}
