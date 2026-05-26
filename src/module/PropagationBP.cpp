#include "crpropa/module/PropagationBP.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {
	void PropagationBP::tryStep(const Y &y, Y &out, Y &error, double h,
			ParticleState particle, double z) const {
		Y outHelp = dY(y.x, y.u, h/2, z, particle);  // 2 steps with h/2
		Y outCompare = dY(outHelp.x, outHelp.u, h/2, z, particle);

		out = dY(y.x, y.u, h, z, particle);  // 1 step with h

		error = errorEstimation(out.x , outCompare.x , h);
	}

	PropagationBP::Y PropagationBP::dY(Vector3d pos, Vector3d dir, double dt, 
		double z, ParticleState current) const {

		// get E and B field at particle position
		Vector3d B = getBFieldAtPosition(pos, z);
		Vector3d E = getEFieldAtPosition(pos, z);

		// if the velocity and the electric field are zero nothing would happen
		if ((E.getR2()==0) && (current.getVelocity().getR2()==0))
			return Y(pos, dir);

		// get some needed candidate features
		Vector3d vel = dir*current.getVelocity().getR();
		double q = current.getCharge();
		// lorentz factor is between 1 and infinity (but never actually infinity)
		double gamma = current.getLorentzFactor();
		double m = gamma*current.getMass();

		// we first do a half leapfrog step
		pos += vel * dt / 2.;

		// if the E field is 0 we only need to do the classical boris push
		// since we are not changing the absolute value of the candidates velocity
		if (E.getR2()==0){
			// Boris help vectors
			Vector3d t = B * q / 2. / m * dt;
			Vector3d s = t * 2. / (1. + t.dot(t));
			Vector3d v_help;

			// Boris push
			v_help = dir + dir.cross(t);
			dir = dir + v_help.cross(s);
		} else {  // otherwise we need to respect the electic field and relativity
			// Boris help vectors
			Vector3d acc = q*E/2./m*dt;  // it is assumed the velocity change is non relativistic
			Vector3d t = B * q / 2. / m * dt;
			Vector3d s = t * 2. / (1. + t.dot(t));

			// differentiate the case for performance improvement:
			// (the relativistic case goes towards the non-relativistic case)
			if (abs(1/gamma - 1) <= 1.e-3){
				Vector3d v_minus = vel + acc;
				Vector3d v_prime = v_minus + v_minus.cross(t);
				v_prime = v_minus + v_prime.cross(s);  // v_prime -> v_plus
				vel = v_prime + acc;  // final velocity
			} else {  // relativistic
				Vector3d v_minus = 1/(1 + acc.dot(vel)/c_squared)*
						( acc/gamma + vel + 1/c_squared*gamma/(1+gamma)*acc.dot(vel)*vel );
				Vector3d v_prime = v_minus + v_minus.cross(t);
				v_prime = v_minus + v_prime.cross(s);  // v_prime -> v_plus
				vel = 1/(1 + acc.dot(v_prime)/c_squared)*
					( acc/gamma + v_prime + 1/c_squared*gamma/(1+gamma)*acc.dot(v_prime)*v_prime );  // final velocity
			}

			// Energy change can only happen when a electric field is present:
			double rm = current.getMass();  // rest mass
			double rm2 = rm*rm;
			double v2 = vel.getR2();
			// dE = E'_kin - E_kin = sqrt(p'^2*c^2 + m^2*c^4) - m*c^2 - E_kin
			deltaE = sqrt(m*m*v2*c_squared + rm2*c_squared*c_squared) - rm*c_squared - current.getEnergy();

			// final velocity might be zero after the the influence of the electric field
			// if that is the case, we know vel must point in the same direction as -acc.
			// since it also has no influence on the next leap frog step we return here:
			if (vel.getR2() == 0)
				return Y(pos, (acc*-1).getUnitVector());
		}

		// the other half leapfrog step (only happens if vel!=0)
		pos += vel * dt / 2.;

		return Y(pos, vel.getUnitVector());
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, double fixedStep) :
		minStep(0) {
		setBField(BField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedStep/c_light);
		setMinimumTimeStep(fixedStep/c_light);
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, double tolerance, double minStep, double maxStep) :
		minStep(0) {
		setBField(BField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxStep/c_light);
		setMinimumTimeStep(minStep/c_light);
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, ref_ptr<ElectricField> EField,
		double fixedStep) : minStep(0)
	{
		setBField(BField);
		setEField(EField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedStep/c_light);
		setMinimumTimeStep(fixedStep/c_light);
	}

	PropagationBP::PropagationBP(ref_ptr<MagneticField> BField, ref_ptr<ElectricField> EField, 
		double tolerance, double minStep,  double maxStep) : minStep(0)
	{
		setBField(BField);
		setEField(EField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxStep/c_light);
		setMinimumTimeStep(minStep/c_light);
	}

	PropagationBP::PropagationBP(double fixedTimeStep, ref_ptr<MagneticField> BField) :
		minStep(0) {
		setBField(BField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedTimeStep);
		setMinimumTimeStep(fixedTimeStep);
	}

	PropagationBP::PropagationBP(double tolerance, double minTimeStep, double maxTimeStep, ref_ptr<MagneticField> BField) :
		minStep(0) {
		setBField(BField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxTimeStep);
		setMinimumTimeStep(minTimeStep);
	}

	PropagationBP::PropagationBP(double fixedTimeStep, ref_ptr<MagneticField> BField, ref_ptr<ElectricField> EField) : 
		minStep(0)
	{
		setBField(BField);
		setEField(EField);
		setTolerance(0.42);
		setMaximumTimeStep(fixedTimeStep);
		setMinimumTimeStep(fixedTimeStep);
	}

	PropagationBP::PropagationBP(double tolerance, double minTimeStep,  double maxTimeStep, 
		ref_ptr<MagneticField> BField, ref_ptr<ElectricField> EField) : minStep(0)
	{
		setBField(BField);
		setEField(EField);
		setTolerance(tolerance);
		setMaximumTimeStep(maxTimeStep);
		setMinimumTimeStep(minTimeStep);
	}

	void PropagationBP::process(Candidate *candidate) const {
		// save the new previous particle state
		ParticleState &current = candidate->current;
		candidate->previous = current;

		Y yIn(current.getPosition(), current.getDirection());

		// calculate charge of particle
		double q = current.getCharge();
		double step = maxStep;

		// rectilinear propagation for neutral particles
		if (q == 0) {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			current.setPosition(yIn.x + yIn.u * candidate->getVelocity() * step);
			candidate->setCurrentStep(step);
			candidate->setNextStep(maxStep);
			return;
		}

		Y yOut, yErr;
		double newStep = step;
		double z = candidate->getRedshift();

		// if minStep is the same as maxStep the adaptive algorithm with its error
		// estimation is not needed and the computation time can be saved:
		if (minStep == maxStep){
			yOut = dY(yIn.x, yIn.u, step, z, current);
		} else {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			newStep = step;
			double r = 42;  // arbitrary value

			// try performing step until the target error (tolerance) or the minimum/maximum step size has been reached
			while (true) {
				tryStep(yIn, yOut, yErr, step, current, z);
				r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
				if (r > 1) {  // large direction error relative to tolerance, try to decrease step size
					if (step == minStep)  // already minimum step size
						break;
					else {
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::max(newStep, 0.1 * step); // limit step size decrease
						newStep = std::max(newStep, minStep); // limit step size to minStep
						step = newStep;
					}
				} else {  // small direction error relative to tolerance, try to increase step size
					if (step != maxStep) {  // only update once if maximum step size yet not reached
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::min(newStep, 5 * step); // limit step size increase
						newStep = std::min(newStep, maxStep); // limit step size to maxStep
					}
					break;
				}
			}
		}

		current.setPosition(yOut.x);
		current.setDirection(yOut.u);
		candidate->setCurrentStep(step);
		candidate->setNextStep(newStep);
		// correct energy, previous already saved by PropagationBP::process
		// candidate->current.setEnergy(candidate->current.getEnergy() + deltaE);
		// deltaE = 0;
	}

	void PropagationBP::setField(ref_ptr<MagneticField> f) {
		setBField(f);
	}

	void PropagationBP::setBField(ref_ptr<MagneticField> BField) {
		this->BField = BField;
	}

	void PropagationBP::setEField(ref_ptr<ElectricField> EField) {
		this->EField = EField;
	}

	Vector3d PropagationBP::getFieldAtPosition(Vector3d pos, double z) const {
		return getBFieldAtPosition(pos, z);
	}

	Vector3d PropagationBP::getBFieldAtPosition(Vector3d pos, double z) const{
		Vector3d B(0, 0, 0);
		try {
			// check if field is valid and use the field vector at the
			// position pos with the redshift z
			if (BField.valid())
				B = BField->getField(pos, z);
		} catch (std::exception &e) {
			KISS_LOG_ERROR 	<< "PropagationBP: Exception in PropagationBP::getFieldAtPosition.\n"
					<< e.what();
		}	
		return B;
	}

	Vector3d PropagationBP::getEFieldAtPosition(Vector3d pos, double z) const{
		Vector3d E(0, 0, 0);
		try {
			// check if field is valid and use the field vector at the
			// position pos with the redshift z
			if (EField.valid())
				E = EField->getField(pos, z);
		} catch (std::exception &e) {
			KISS_LOG_ERROR 	<< "PropagationBP: Exception in PropagationBP::getFieldAtPosition.\n"
					<< e.what();
		}	
		return E;
	}

	double PropagationBP::errorEstimation(const Vector3d x1, const Vector3d x2, double step) const {
		// compare the position after one step with the position after two steps with step/2.
		Vector3d diff = (x1 - x2);

		double S = diff.getR() / (step * (1 - 1/4.) );	// 1/4 = (1/2)²  number of steps for x1 divided by number of steps for x2 to the power of p (order)

		return S;
	}

	void PropagationBP::setTolerance(double tol) {
		if ((tol > 1) or (tol < 0))
			throw std::runtime_error(
					"PropagationBP: target error not in range 0-1");
		tolerance = tol;
	}

	void PropagationBP::setMinimumStep(double min) {
		if (min/c_light < 0)
			throw std::runtime_error("PropagationBP: minStep < 0 ");
		if (min/c_light > maxStep)
			throw std::runtime_error("PropagationBP: minStep > maxStep");
		minStep = min/c_light;
	}

	void PropagationBP::setMaximumStep(double max) {
		if (max/c_light < minStep)
			throw std::runtime_error("PropagationBP: maxStep < minStep");
		maxStep = max/c_light;
	}

	void PropagationBP::setMinimumTimeStep(double min) {
		if (min < 0)
			throw std::runtime_error("PropagationBP: minStep < 0 ");
		if (min > maxStep)
			throw std::runtime_error("PropagationBP: minStep > maxStep");
		minStep = min;
	}

	void PropagationBP::setMaximumTimeStep(double max) {
		if (max < minStep)
			throw std::runtime_error("PropagationBP: maxStep < minStep");
		maxStep = max;
	}

	std::string PropagationBP::getDescription() const {
		std::stringstream s;
		s << "Propagation in magnetic fields using the adaptive Boris push method.";
		s << " Target error: " << tolerance;
		s << ", Minimum Step: " << minStep / kiloyear << " kiloyear";
		s << ", Maximum Step: " << maxStep / kiloyear << " kiloyear";
		return s.str();
	}
} // namespace crpropa
