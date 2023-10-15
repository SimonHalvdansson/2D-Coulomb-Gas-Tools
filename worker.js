
self.addEventListener('message', function(e) {
  const i = e.data.index;
  const particles = new Float64Array(e.data.bufferToSend);

  const n = particles.length/4;
  const h = e.data.h;
  const dt = e.data.dt;
  const pot = e.data.pot;

  function Q(x, y) {
    q = 0

    //Ginibre
    if (pot == 0)
      q += x*x+y*y
    //Mittag-Leffler
    if (pot == 1)
      q += Math.pow(x*x+y*y, 2)
    //mittag-leffler lambda=10
    if (pot == 2)
      q += Math.pow(x*x+y*y, 10)

    return q
    
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