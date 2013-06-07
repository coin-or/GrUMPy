# BAK_visual.py
#
#  Copyright 2009 Google Inc.
#  Copyright 2007 University of Pittsburgh.
#  Google coding done by Brady Hunsaker.
#  U of Pittsburgh coding done by Osman Ozaltin and Brady Hunsaker.
#
#  This file is part of BAK (Branch-and-bound Analysis Kit).
#
#  The contents of this file are subject to the Common Public License
#  1.0.  (the "License"); you may not use this file except in
#  compliance with the License. You should have received a copy of
#  the Common Public License along with STOP.
#
#  Software distributed under the License is distributed on an "AS
#  IS" basis, WITHOUT WARRANTY OF ANY KIND, either express or
#  implied. See the License for the specific language governing
#  rights and limitations under the License.
#
#  Alternatively, the contents of this file may be used under the
#  terms of the GNU General Public License Version 2 or later (the
#  "GPL"), in which case the provisions of the GPL are applicable
#  instead of those above. If you wish to allow use of your version
#  of this file only under the terms of the GPL, and not to allow
#  others to use your version of this file under the terms of the
#  CPL, indicate your decision by deleting the provisions above and
#  replace them with the notice and other provisions required by the
#  GPL. If you do not delete the provisions above, a recipient may
#  use your version of this file under the terms of either the CPL or
#  the GPL.
#

# For developers: Please keep the code style consistent with the Python
# style guide for Google's Summer of Code, except use 4 spaces to indent:
#   http://code.google.com/p/soc/wiki/PythonStyleGuide

__author__ = 'Brady Hunsaker, Osman Ozaltin'
__maintainer__ = 'Brady Hunsaker (bhunsaker@google.com)'

"""Double exponential smoothing forecasting in chained sequences.

A single sequence is stored in a DoubleExponentialSmoothingForecaster, which
also provides forecasts for time to completion based on the stored measures,
which should be monotonically decreasing.

Sequences are chained together in a ForecastingChainedSequences object.
Such an object includes scale factors for each sequence.  The scale factors
will be applied to the measurements, usually for the purpose of making the
chained sequence monotonically decreasing.
"""

import math


class ProgressMeasurement(object):
    """Key data recording progress.  Data members are public.
    """
    def __init__(self, time, value, active_node_count, node_count):
        self.time = float(time)
        self.value = float(value)
        self.active_node_count = int(active_node_count)
        self.node_count = int(node_count)


class TimeForecast(object):
    """Time-stamped forecast of time remaining.
    """
    def __init__(self, time, forecast):
        self.time = float(time)
        self.forecast = float(forecast)


class DoubleExponentialSmoothingForecaster(object):
    """Uses double exponential smoothing to forecast values.
    """
    def __init__(self, scale_factor, first_value, first_time):
        self._measures = []
        self._forecasts = []

        self._b_t = None
        self._S_t = None

        self._scale_factor = float(scale_factor)
        self._first_time = -1

        if first_value is None:
            self._first_value = 1
        else:
            self._first_value = float(first_value)

        self._alpha = 0.5
        self._lambda = 0.5
        self._gamma = 0.5

        self._counter = 0

    def AddMeasure(self, measurement):
        self._measures.append(measurement)
        self.ComputeForecast()

    def ComputeForecast(self):
        # Check if the forecast was already computed
        if len(self._forecasts) == len(self._measures) - 3:
            return
        # Must have at least 4 measurements to make a forecast.
        if len(self._measures) < 4:
            return

        # Check that a minimum amount of time has passed
        if self._measures[-1].time - self._measures[-2].time < 1:
            return

        measure_value = self._measures[-1].value * self._scale_factor
        previous_value = self._measures[-2].value * self._scale_factor
        time = self._measures[-1].time
        previous_time = self._measures[-2].time
        delta = float((measure_value - previous_value)) / float((time - previous_time))
        if -delta < 0.00001 * measure_value:
            delta = 0

        # Set initial parameters if necessary otherwise update the estimators
        if self._b_t is None:
            if self._first_value > 0:
                self._b_t = delta
                self._S_t = (measure_value * self._alpha +
                             previous_value * (1 - self._alpha))
            else:
                self._b_t = self._first_value
                self._S_t = (measure_value * self._alpha +
                             previous_value * (1 - self._alpha))

        if self._b_t < -measure_value:
            self._b_t = -measure_value
            self._S_t = (measure_value * self._alpha +
                         previous_value * (1 - self._alpha))

        updated = False

        if delta < 0:
            self._b_t = self._gamma * delta + (1 - self._gamma) * self._b_t
            self._S_t = (self._alpha * measure_value +
                         (1 - self._alpha) * self._S_t)
            updated = True

        print 'delta: %f' %delta
        print 'self._lambda*self._b_t: %f' %(self._lambda*self._b_t)
        if delta < self._lambda * self._b_t:
            if not updated:
                self._b_t = self._gamma * delta + (1 - self._gamma) * self._b_t
                self._S_t = (self._alpha * measure_value +
                             (1 - self._alpha) * self._S_t)
                updated = True

            forecast = (time + (float(-self._S_t) / min(delta,self._b_t)))
            print 'A', time, previous_value, measure_value,
            print self._S_t, self._b_t, forecast
            #cent = float(input("STOP"))
        elif len(self._forecasts) >= 1:
            # The measure didn't change but we have a previous forecast
            if self._forecasts[-1].forecast >= time:
                forecast = self._forecasts[-1].forecast
            else:
                forecast = (time +
                            (float(self._measures[-1].active_node_count) /
                             self._measures[-2].active_node_count) *
                            (self._forecasts[-1].forecast -
                             self._forecasts[-1].time))

            print 'B', time, previous_value, measure_value, forecast
            #cent = float(input("STOP"))
        else:
            # The measure didn't change and we have no previous forecast
            forecast = (time +
                        float(self._measures[-1].active_node_count * time) /
                        (self._measures[-1].node_count -
                         self._measures[-1].active_node_count))
            print 'C', time, previous_value, measure_value, forecast

        self._forecasts.append(TimeForecast(time, forecast))

        #if (time >= 200):
        #    cent = float(input("STOP"))

        #if (forecast >= 60000):
        #    cent = float(input("STOP"))


    def GetForecasts(self):
        return self._forecasts

    def GetMeasures(self):
        return self._measures


class ForecastingChainedSequences(object):
    """
    """
    def __init__(self):
        self._sequences = []
        self._scale_factors = []
        self._current_sequence = -1

    def AddMeasure(self, time, value, active_node_count, node_count):
        assert self._current_sequence >= 0
        self._sequences[self._current_sequence].AddMeasure(
            ProgressMeasurement(time, value, active_node_count, node_count))

    def StartNewSequence(self, scale_factor):
        """Starts a new sequence of measures.

        The scale factor should compare only to the previous sequence.
        """
        if self._scale_factors:
            scale_factor *= self._scale_factors[-1]
        self._scale_factors.append(scale_factor)

        assert len(self._sequences) == self._current_sequence + 1

        if self._current_sequence < 0:
            self._sequences.append(DoubleExponentialSmoothingForecaster(
                scale_factor, -1, -1))
        else:
            self._sequences.append(DoubleExponentialSmoothingForecaster(
                scale_factor, self._sequences[self._current_sequence]._b_t,
                self._sequences[self._current_sequence]._first_time))

        self._current_sequence += 1

        assert len(self._scale_factors) == len(self._sequences)

    def GetAllForecasts(self):
        forecasts = []
        for sequence in self._sequences:
            forecasts.extend(sequence.GetForecasts())
        return forecasts

    def GetAllMeasures(self):
        measures = []
        for i, sequence in enumerate(self._sequences):
            temp_measures = sequence.GetMeasures()
            for measure in temp_measures:
                measure.value *= self._scale_factors[i]
            measures.extend(temp_measures)
        return measures
