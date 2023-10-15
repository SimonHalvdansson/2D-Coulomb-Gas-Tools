
self.addEventListener('message', function(e) {
  const i = e.data.index;
  const particles = new Float64Array(e.data.bufferToSend);

  const n = particles.length/4;
  const h = e.data.h;
  const dt = e.data.dt;
  const pot = e.data.pot;
  const lemniscateT = e.data.lemniscateT;
  const lemInterpol = e.data.lemInterpol;
  const ca = e.data.ca;
  const cb = e.data.cb;
  const bx = e.data.bx;
  const by = e.data.by;
  const ax = e.data.ax;
  const ay = e.data.ay;

  function Q(x, y) {
    let q = 0;
    let real = x;
    let imag = y;

    let r = Math.sqrt(real * real + imag * imag);
    let theta = Math.atan2(imag, real);

    if (pot === 0)
      q += x * x + y * y;
    if (pot === 1)
      q += Math.pow(x * x + y * y, 2);
    if (pot === 2)
      q += Math.pow(x * x + y * y, 10);
    if (pot === 3)
      q += Math.pow(x * x + y * y, 2) - lemniscateT * 2 / Math.sqrt(2) * Math.pow(r, 2) * Math.cos(2 * theta);
    if (pot === 4)
      q += Math.pow(x * x + y * y, 3) - lemniscateT * 2 / Math.sqrt(3) * Math.pow(r, 3) * Math.cos(3 * theta);
    if (pot === 5)
      q += Math.pow(x * x + y * y, 3) - lemniscateT * 2 / Math.sqrt(5) * Math.pow(r, 5) * Math.cos(5 * theta);
    if (pot === 6)
      q += Math.pow(x * x + y * y, lemInterpol) - lemniscateT * 2 / Math.sqrt(lemInterpol) * Math.pow(r, lemInterpol) * Math.cos(lemInterpol * theta);

    q += (-cb * Math.log(Math.sqrt(Math.pow(x - bx, 2) + Math.pow(y - by, 2))) - ca * Math.log(Math.sqrt(Math.pow(x - ax, 2) + Math.pow(y + ay, 2))));

    return q;
  }

  //first do the Q derivative numerically (as it is not O(n^2))
  let acc_x = -n*(Q(particles[i*4+0]+h, particles[i*4+1]) - Q(particles[i*4+0], particles[i*4+1]))/h
  let acc_y = -n*(Q(particles[i*4+0], particles[i*4+1]+h) - Q(particles[i*4+0], particles[i*4+1]))/h

  /*
  now for the pairwise interaction, the force should be in the
  direction of the difference vector between the two particles
  and be proportional to the derivative of the energy
  */ 

  for (let j = 0; j < n; j++) {
    if (j != i) {

      let px = particles[i*4+0]
      let py = particles[i*4+1]
      let qx = particles[j*4+0]
      let qy = particles[j*4+1]

      //for now, let's just do numerical derivative
      let den = (px*px - 2*px*qx + qx*qx + py*py - 2*py*qy + qy*qy);
      acc_x -= (2*(qx - px))/den;
      acc_y -= (2*(qy - py))/den;
    }
  }

  self.postMessage({ i, acc_x, acc_y });
});