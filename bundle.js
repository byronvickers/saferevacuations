(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){

lineIntersect = require('@turf/line-intersect').default
lineDistance = require('@turf/point-to-line-distance').default
distance = require('@turf/distance').default
centroid = require('@turf/centroid').default
agentmaps = require('agentmaps').default
routing = require('agentmaps/src/routing')

////// dario code ////

map_data.features = map_data.features.filter((item) => item.properties.tags.highway && ['trunk', 'primary', 'secondary', 'tertiary', 'unclassified'].indexOf(item.properties.tags.highway) !== -1)
map = L.map("demo_map").fitBounds(bounding_points);
map.setZoom(map.getZoom()+1)
L.tileLayer(
    "http://{s}.tiles.wmflabs.org/bw-mapnik/{z}/{x}/{y}.png",
    // "http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
    {
        attribution: "Thanks to <a href=\"http://openstreetmap.org\">OpenStreetMap</a> community",
    }
).addTo(map);

let street_options = {
  "color": "yellow",
  "weight": 4,
  "opacity": .5
};
let unit_options = {
  front_buffer: 6,
  side_buffer: 3,
  length: 10,
  depth: 18,
  color: "blue"
};
// let embassies_coords = embassies.map((item) => {
//   return [item.geometry.y, item.geometry.x]
// });
// embassies_coords = embassies_coords.filter((item) => item[0] > bounding_points[0][0] && item[0] < bounding_points[1][0] && item[1] > bounding_points[1][1] && item[1] < bounding_points[1][1])


agentmap = L.A.agentmap(map);

let options_set = {
  school: {
    "color": "#e34a33",
    "weight": 1,
    "opacity": .87,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 100,
    "depth": 100
  },
  embassy: {
    "color": "green",
    "weight": 1,
    "opacity": .87,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 100,
    "depth": 100
  },
  standard: {
    "color": "grey",
    "weight": 1,
    "opacity": .3,
    "front_buffer": 6,
    "side_buffer": 3,
    "length": 25,
    "depth": 25
  },
  person: {
    radius: .5,
    color: "red",
    fillColor: "red",
    opacity: 0.7
  },
  student: {
    radius: .5,
    color: "green",
    fillColor: "green",
    opacity: 0.7
  },
  embassyVip: {
    radius: .5,
    color: "blue",
    fillColor: "blue",
    opacity: 0.7
  }
};

unit_options = options_set.standard
congregate_idxs = [0,1,2]

agentmap.buildingify(bounding_points, map_data, street_options, unit_options, units_data);
congregate_idxs.forEach(function(idx){
    unit = agentmap.units.getLayers()[idx]
    unit.feature.properties.isCongregationPoint = true
    circ = new L.circle(unit.getCenter(), 1, {
      fillOpacity: 0.0,
      opacity: 0.0
    })
    .bindTooltip('EVAC', { permanent: true, direction: "center", offset: [0,10], className: "flagstyle"})
    .addTo(map);
  }
)
agentmap.units.eachLayer(function(unit){
  if(false){//unit.feature.properties.isCongregationPoint){
    //
  } else {
    unittype = unit.feature.properties.poi_type
    unitmarker = false
    if(unittype === "embassy"){
      unitmarker = '🤝🏾'
    } else if (unittype === "school") {
      unitmarker = '🎓'
    }
    if(unitmarker){
      circ = new L.circle(unit.getCenter(), 1, {
        fillOpacity: 0.0,
        opacity: 0.0
      })
      .bindTooltip(unitmarker, { permanent: true, direction: "center", offset: [-10,0], className: "labelstyle"})
      .addTo(map);
    }
  }
})



/////

agentmap.units.eachLayer(function(unit){
  unit.options = Object.assign(unit.options, unit.feature.options)
  // unit.options.fillColor = unit.feature.options
  unit.setStyle()
})

////////

function destroygrid(grid) {
  if (map.hasLayer(grid)) {
    grid.unregister();
    map.removeLayer(grid);
  }
}

/* seeded js rng from https://stackoverflow.com/questions/521295/seeding-the-random-number-generator-in-javascript */
var seed = 1;
function random() {
    var x = Math.sin(seed++) * 10000;
    return x - Math.floor(x);
}

/* Durstenfeld shuffle from https://stackoverflow.com/questions/2450954/how-to-randomize-shuffle-a-javascript-array */
function shuffleArray(array) {
    for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(random() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
    }
}

function randomUnitAgentMaker(id){
	let index = this.agents.count();

  sel = Math.random()
  if(sel < 0.2){
    poitype = "embassy"
  } else if (sel < 0.4) {
    poitype = "school"
  } else {
    poitype = "standard"
  }

  unitset = this.units.getLayers().filter(unit => unit.feature.properties.poi_type === poitype)

  // let unit = this.units.getLayers()[Math.floor(agentmap.units.count() * Math.random())],
  let unit = unitset[Math.floor(unitset.length * Math.random())],
	unit_id = this.units.getLayerId(unit),
	center_point = centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};

  if (unit.feature.properties.poi_type === 'school') {
    center_point.properties.layer_options = Math.random() > 0.5 ? options_set.person : options_set.student;
  } else if (unit.feature.properties.poi_type === 'embassy') {
    center_point.properties.layer_options = Math.random() > 0.8 ? options_set.person : options_set.embassyVip;
  } else {
    center_point.properties.layer_options = options_set.person
  }

  return center_point;

}


window.numAgents = 500
idxarr = [...Array(numAgents).keys()];
shuffleArray(idxarr)
agentmap.agentify(numAgents, randomUnitAgentMaker);
// agentmap.agentify(numAgents, agentmap.seqUnitAgentMaker);

window.scatter = false
window.commandwaiting = false

agentmap.controller = function() {
    if ((agentmap.state.ticks === 0) || (window.commandwaiting)) {
        agentmap.agents.eachLayer(function(agent) {
              if(window.scatter){
                dest_index = Math.floor(agentmap.units.count() * Math.random())
                dest_unit = agentmap.units.getLayers()[dest_index]
              } else {
                closest_unit = null
                closest_dist = 9999999
                agentgeo = agent.toGeoJSON()
                congregate_idxs.forEach(function(idx){
                  unit = agentmap.units.getLayers()[idx]
                  unitgeo = centroid(unit.feature)
                  if(distance(agentgeo, unitgeo) < closest_dist){
                    closest_unit = unit
                    closest_dist = distance(agentgeo, unitgeo)
                  }
                })
                dest_unit = closest_unit
              }

              dest_unit_id = agentmap.units.getLayerId(dest_unit)
              dest_unit_center = dest_unit.getBounds().getCenter()

              if(agent.trip.path.length > 0){
                closest_unit = null
                closest_dist = 9999999
                agentgeo = agent.toGeoJSON()
                agentmap.units.eachLayer(function(unit){
                  unitgeo = centroid(unit.feature)
                  if(distance(agentgeo, unitgeo) < closest_dist){
                    closest_unit = unit
                    closest_dist = distance(agentgeo, unitgeo)
                  }
                })
                closest_unit_center = closest_unit.getBounds().getCenter()
                agent.setTravelToPlace(closest_unit_center, {type: "unit", id: agentmap.units.getLayerId(closest_unit)}, 4, true, true);
              }

              rndspeed = 4*(1 + random()/4)
              rnddelay = Math.floor(Math.random()*10000)
              // sendAgentWithDelay(agentmap.)
              try{
                agent.setTravelToPlace(dest_unit_center, {type: "unit", id: dest_unit_id}, rndspeed)
              }
              catch(error){
                // just leave them be, for now.
              }
        })
        window.commandwaiting = false
    }
    agentmap.agents.eachLayer(function(agent){
      agent.moveIt();
    })
    if (agentmap.state.ticks % 25 === 0){
      if(typeof grid !== 'undefined'){
        destroygrid(grid)
        grid = creategrid()
      }

    }
}
agentmap.run();
agentmap.setAnimationInterval(25)


// agentmap.streets.eachLayer(function(unit){unit.once("click", function(street){agentmap.streets.removeLayer(street)})})
// agentmap.streets.eachLayer(function(street){street.on("click", handleClick)})
map.on('click', function(e){
    clickCircle = new L.circle(e.latlng, 140*Math.pow(2, 15-map.getZoom()), {
      color: '#f07300',
      fillOpacity: 0.3,
      opacity: 0.5
    })
    .bindTooltip('🔥', { permanent: true, direction: "center", className: "labelstyle"})
    .addTo(map);

    // clickCircle.bringToFront()
    // window.c = clickCircle

    street_dists = []
    agentmap.streets.eachLayer(function(street){
      if(lineDistance(clickCircle.toGeoJSON(), street.feature)*1000 < clickCircle.getRadius()*0.8){ /* to get metres */
        agentmap.streets.removeLayer(street)
        street_id = agentmap.streets.getLayerId(street)

        agentmap.units.eachLayer(function(unit){
          if(unit.street_id === street_id){
            agentmap.units.removeLayer(unit)
          }
        });
      }
    })

    agentmap.streets.graph = routing.streetsToGraph(agentmap.streets),
  	agentmap.pathfinder = routing.getPathFinder(agentmap.streets.graph);

    window.commandwaiting = true
});

window.disperse = function(){
  window.scatter = true
  window.commandwaiting = true
}

window.congregate = function(){
  window.scatter = false
  window.commandwaiting = true
}

window.pausetoggle = function(){
  if(agentmap.state.paused){
    agentmap.run()
  } else{
    agentmap.pause()
  }
}

},{"@turf/centroid":10,"@turf/distance":14,"@turf/line-intersect":20,"@turf/point-to-line-distance":38,"agentmaps":46,"agentmaps/src/routing":47}],2:[function(require,module,exports){
"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
}
Object.defineProperty(exports, "__esModule", { value: true });
var bearing_1 = __importDefault(require("@turf/bearing"));
var destination_1 = __importDefault(require("@turf/destination"));
var distance_1 = __importDefault(require("@turf/distance"));
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
/**
 * Takes a {@link LineString} and returns a {@link Point} at a specified distance along the line.
 *
 * @name along
 * @param {Feature<LineString>} line input line
 * @param {number} distance distance along the line
 * @param {Object} [options] Optional parameters
 * @param {string} [options.units="kilometers"] can be degrees, radians, miles, or kilometers
 * @returns {Feature<Point>} Point `distance` `units` along the line
 * @example
 * var line = turf.lineString([[-83, 30], [-84, 36], [-78, 41]]);
 * var options = {units: 'miles'};
 *
 * var along = turf.along(line, 200, options);
 *
 * //addToMap
 * var addToMap = [along, line]
 */
function along(line, distance, options) {
    if (options === void 0) { options = {}; }
    // Get Coords
    var geom = invariant_1.getGeom(line);
    var coords = geom.coordinates;
    var travelled = 0;
    for (var i = 0; i < coords.length; i++) {
        if (distance >= travelled && i === coords.length - 1) {
            break;
        }
        else if (travelled >= distance) {
            var overshot = distance - travelled;
            if (!overshot) {
                return helpers_1.point(coords[i]);
            }
            else {
                var direction = bearing_1.default(coords[i], coords[i - 1]) - 180;
                var interpolated = destination_1.default(coords[i], overshot, direction, options);
                return interpolated;
            }
        }
        else {
            travelled += distance_1.default(coords[i], coords[i + 1], options);
        }
    }
    return helpers_1.point(coords[coords.length - 1]);
}
exports.default = along;

},{"@turf/bearing":4,"@turf/destination":13,"@turf/distance":14,"@turf/helpers":15,"@turf/invariant":17}],3:[function(require,module,exports){
'use strict';

var meta = require('@turf/meta');

/**
 * Takes a set of features, calculates the bbox of all input features, and returns a bounding box.
 *
 * @name bbox
 * @param {GeoJSON} geojson any GeoJSON object
 * @returns {BBox} bbox extent in [minX, minY, maxX, maxY] order
 * @example
 * var line = turf.lineString([[-74, 40], [-78, 42], [-82, 35]]);
 * var bbox = turf.bbox(line);
 * var bboxPolygon = turf.bboxPolygon(bbox);
 *
 * //addToMap
 * var addToMap = [line, bboxPolygon]
 */
function bbox(geojson) {
    var BBox = [Infinity, Infinity, -Infinity, -Infinity];
    meta.coordEach(geojson, function (coord) {
        if (BBox[0] > coord[0]) BBox[0] = coord[0];
        if (BBox[1] > coord[1]) BBox[1] = coord[1];
        if (BBox[2] < coord[0]) BBox[2] = coord[0];
        if (BBox[3] < coord[1]) BBox[3] = coord[1];
    });
    return BBox;
}

module.exports = bbox;
module.exports.default = bbox;

},{"@turf/meta":34}],4:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
// http://en.wikipedia.org/wiki/Haversine_formula
// http://www.movable-type.co.uk/scripts/latlong.html
/**
 * Takes two {@link Point|points} and finds the geographic bearing between them,
 * i.e. the angle measured in degrees from the north line (0 degrees)
 *
 * @name bearing
 * @param {Coord} start starting Point
 * @param {Coord} end ending Point
 * @param {Object} [options={}] Optional parameters
 * @param {boolean} [options.final=false] calculates the final bearing if true
 * @returns {number} bearing in decimal degrees, between -180 and 180 degrees (positive clockwise)
 * @example
 * var point1 = turf.point([-75.343, 39.984]);
 * var point2 = turf.point([-75.534, 39.123]);
 *
 * var bearing = turf.bearing(point1, point2);
 *
 * //addToMap
 * var addToMap = [point1, point2]
 * point1.properties['marker-color'] = '#f00'
 * point2.properties['marker-color'] = '#0f0'
 * point1.properties.bearing = bearing
 */
function bearing(start, end, options) {
    if (options === void 0) { options = {}; }
    // Reverse calculation
    if (options.final === true) {
        return calculateFinalBearing(start, end);
    }
    var coordinates1 = invariant_1.getCoord(start);
    var coordinates2 = invariant_1.getCoord(end);
    var lon1 = helpers_1.degreesToRadians(coordinates1[0]);
    var lon2 = helpers_1.degreesToRadians(coordinates2[0]);
    var lat1 = helpers_1.degreesToRadians(coordinates1[1]);
    var lat2 = helpers_1.degreesToRadians(coordinates2[1]);
    var a = Math.sin(lon2 - lon1) * Math.cos(lat2);
    var b = Math.cos(lat1) * Math.sin(lat2) -
        Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);
    return helpers_1.radiansToDegrees(Math.atan2(a, b));
}
/**
 * Calculates Final Bearing
 *
 * @private
 * @param {Coord} start starting Point
 * @param {Coord} end ending Point
 * @returns {number} bearing
 */
function calculateFinalBearing(start, end) {
    // Swap start & end
    var bear = bearing(end, start);
    bear = (bear + 180) % 360;
    return bear;
}
exports.default = bearing;

},{"@turf/helpers":15,"@turf/invariant":17}],5:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var invariant_1 = require("@turf/invariant");
// http://en.wikipedia.org/wiki/Even%E2%80%93odd_rule
// modified from: https://github.com/substack/point-in-polygon/blob/master/index.js
// which was modified from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
/**
 * Takes a {@link Point} and a {@link Polygon} or {@link MultiPolygon} and determines if the point
 * resides inside the polygon. The polygon can be convex or concave. The function accounts for holes.
 *
 * @name booleanPointInPolygon
 * @param {Coord} point input point
 * @param {Feature<Polygon|MultiPolygon>} polygon input polygon or multipolygon
 * @param {Object} [options={}] Optional parameters
 * @param {boolean} [options.ignoreBoundary=false] True if polygon boundary should be ignored when determining if
 * the point is inside the polygon otherwise false.
 * @returns {boolean} `true` if the Point is inside the Polygon; `false` if the Point is not inside the Polygon
 * @example
 * var pt = turf.point([-77, 44]);
 * var poly = turf.polygon([[
 *   [-81, 41],
 *   [-81, 47],
 *   [-72, 47],
 *   [-72, 41],
 *   [-81, 41]
 * ]]);
 *
 * turf.booleanPointInPolygon(pt, poly);
 * //= true
 */
function booleanPointInPolygon(point, polygon, options) {
    if (options === void 0) { options = {}; }
    // validation
    if (!point) {
        throw new Error("point is required");
    }
    if (!polygon) {
        throw new Error("polygon is required");
    }
    var pt = invariant_1.getCoord(point);
    var geom = invariant_1.getGeom(polygon);
    var type = geom.type;
    var bbox = polygon.bbox;
    var polys = geom.coordinates;
    // Quick elimination if point is not inside bbox
    if (bbox && inBBox(pt, bbox) === false) {
        return false;
    }
    // normalize to multipolygon
    if (type === "Polygon") {
        polys = [polys];
    }
    var insidePoly = false;
    for (var i = 0; i < polys.length && !insidePoly; i++) {
        // check if it is in the outer ring first
        if (inRing(pt, polys[i][0], options.ignoreBoundary)) {
            var inHole = false;
            var k = 1;
            // check for the point in any of the holes
            while (k < polys[i].length && !inHole) {
                if (inRing(pt, polys[i][k], !options.ignoreBoundary)) {
                    inHole = true;
                }
                k++;
            }
            if (!inHole) {
                insidePoly = true;
            }
        }
    }
    return insidePoly;
}
exports.default = booleanPointInPolygon;
/**
 * inRing
 *
 * @private
 * @param {Array<number>} pt [x,y]
 * @param {Array<Array<number>>} ring [[x,y], [x,y],..]
 * @param {boolean} ignoreBoundary ignoreBoundary
 * @returns {boolean} inRing
 */
function inRing(pt, ring, ignoreBoundary) {
    var isInside = false;
    if (ring[0][0] === ring[ring.length - 1][0] && ring[0][1] === ring[ring.length - 1][1]) {
        ring = ring.slice(0, ring.length - 1);
    }
    for (var i = 0, j = ring.length - 1; i < ring.length; j = i++) {
        var xi = ring[i][0];
        var yi = ring[i][1];
        var xj = ring[j][0];
        var yj = ring[j][1];
        var onBoundary = (pt[1] * (xi - xj) + yi * (xj - pt[0]) + yj * (pt[0] - xi) === 0) &&
            ((xi - pt[0]) * (xj - pt[0]) <= 0) && ((yi - pt[1]) * (yj - pt[1]) <= 0);
        if (onBoundary) {
            return !ignoreBoundary;
        }
        var intersect = ((yi > pt[1]) !== (yj > pt[1])) &&
            (pt[0] < (xj - xi) * (pt[1] - yi) / (yj - yi) + xi);
        if (intersect) {
            isInside = !isInside;
        }
    }
    return isInside;
}
/**
 * inBBox
 *
 * @private
 * @param {Position} pt point [x,y]
 * @param {BBox} bbox BBox [west, south, east, north]
 * @returns {boolean} true/false if point is inside BBox
 */
function inBBox(pt, bbox) {
    return bbox[0] <= pt[0] &&
        bbox[1] <= pt[1] &&
        bbox[2] >= pt[0] &&
        bbox[3] >= pt[1];
}

},{"@turf/invariant":17}],6:[function(require,module,exports){
'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var center = _interopDefault(require('@turf/center'));
var turfBbox = _interopDefault(require('@turf/bbox'));
var turfJsts = require('turf-jsts');
var projection = require('@turf/projection');
var meta = require('@turf/meta');
var d3Geo = require('d3-geo');
var helpers = require('@turf/helpers');

/**
 * Calculates a buffer for input features for a given radius. Units supported are miles, kilometers, and degrees.
 *
 * When using a negative radius, the resulting geometry may be invalid if
 * it's too small compared to the radius magnitude. If the input is a
 * FeatureCollection, only valid members will be returned in the output
 * FeatureCollection - i.e., the output collection may have fewer members than
 * the input, or even be empty.
 *
 * @name buffer
 * @param {FeatureCollection|Geometry|Feature<any>} geojson input to be buffered
 * @param {number} radius distance to draw the buffer (negative values are allowed)
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units="kilometers"] any of the options supported by turf units
 * @param {number} [options.steps=64] number of steps
 * @returns {FeatureCollection|Feature<Polygon|MultiPolygon>|undefined} buffered features
 * @example
 * var point = turf.point([-90.548630, 14.616599]);
 * var buffered = turf.buffer(point, 500, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [point, buffered]
 */
function buffer(geojson, radius, options) {
    // Optional params
    options = options || {};
    var units = options.units;
    var steps = options.steps || 64;

    // validation
    if (!geojson) throw new Error('geojson is required');
    if (typeof options !== 'object') throw new Error('options must be an object');
    if (typeof steps !== 'number') throw new Error('steps must be an number');

    // Allow negative buffers ("erosion") or zero-sized buffers ("repair geometry")
    if (radius === undefined) throw new Error('radius is required');
    if (steps <= 0) throw new Error('steps must be greater than 0');

    // default params
    steps = steps || 64;
    units = units || 'kilometers';

    var results = [];
    switch (geojson.type) {
    case 'GeometryCollection':
        meta.geomEach(geojson, function (geometry) {
            var buffered = bufferFeature(geometry, radius, units, steps);
            if (buffered) results.push(buffered);
        });
        return helpers.featureCollection(results);
    case 'FeatureCollection':
        meta.featureEach(geojson, function (feature$$1) {
            var multiBuffered = bufferFeature(feature$$1, radius, units, steps);
            if (multiBuffered) {
                meta.featureEach(multiBuffered, function (buffered) {
                    if (buffered) results.push(buffered);
                });
            }
        });
        return helpers.featureCollection(results);
    }
    return bufferFeature(geojson, radius, units, steps);
}

/**
 * Buffer single Feature/Geometry
 *
 * @private
 * @param {Feature<any>} geojson input to be buffered
 * @param {number} radius distance to draw the buffer
 * @param {string} [units='kilometers'] any of the options supported by turf units
 * @param {number} [steps=64] number of steps
 * @returns {Feature<Polygon|MultiPolygon>} buffered feature
 */
function bufferFeature(geojson, radius, units, steps) {
    var properties = geojson.properties || {};
    var geometry = (geojson.type === 'Feature') ? geojson.geometry : geojson;

    // Geometry Types faster than jsts
    if (geometry.type === 'GeometryCollection') {
        var results = [];
        meta.geomEach(geojson, function (geometry) {
            var buffered = bufferFeature(geometry, radius, units, steps);
            if (buffered) results.push(buffered);
        });
        return helpers.featureCollection(results);
    }

    // Project GeoJSON to Transverse Mercator projection (convert to Meters)
    var projected;
    var bbox = turfBbox(geojson);
    var needsTransverseMercator = bbox[1] > 50 && bbox[3] > 50;

    if (needsTransverseMercator) {
        projected = {
            type: geometry.type,
            coordinates: projectCoords(geometry.coordinates, defineProjection(geometry))
        };
    } else {
        projected = projection.toMercator(geometry);
    }

    // JSTS buffer operation
    var reader = new turfJsts.GeoJSONReader();
    var geom = reader.read(projected);
    var distance = helpers.radiansToLength(helpers.lengthToRadians(radius, units), 'meters');
    var buffered = turfJsts.BufferOp.bufferOp(geom, distance);
    var writer = new turfJsts.GeoJSONWriter();
    buffered = writer.write(buffered);

    // Detect if empty geometries
    if (coordsIsNaN(buffered.coordinates)) return undefined;

    // Unproject coordinates (convert to Degrees)
    var result;
    if (needsTransverseMercator) {
        result = {
            type: buffered.type,
            coordinates: unprojectCoords(buffered.coordinates, defineProjection(geometry))
        };
    } else {
        result = projection.toWgs84(buffered);
    }

    return (result.geometry) ? result : helpers.feature(result, properties);
}

/**
 * Coordinates isNaN
 *
 * @private
 * @param {Array<any>} coords GeoJSON Coordinates
 * @returns {boolean} if NaN exists
 */
function coordsIsNaN(coords) {
    if (Array.isArray(coords[0])) return coordsIsNaN(coords[0]);
    return isNaN(coords[0]);
}

/**
 * Project coordinates to projection
 *
 * @private
 * @param {Array<any>} coords to project
 * @param {GeoProjection} proj D3 Geo Projection
 * @returns {Array<any>} projected coordinates
 */
function projectCoords(coords, proj) {
    if (typeof coords[0] !== 'object') return proj(coords);
    return coords.map(function (coord) {
        return projectCoords(coord, proj);
    });
}

/**
 * Un-Project coordinates to projection
 *
 * @private
 * @param {Array<any>} coords to un-project
 * @param {GeoProjection} proj D3 Geo Projection
 * @returns {Array<any>} un-projected coordinates
 */
function unprojectCoords(coords, proj) {
    if (typeof coords[0] !== 'object') return proj.invert(coords);
    return coords.map(function (coord) {
        return unprojectCoords(coord, proj);
    });
}

/**
 * Define Transverse Mercator projection
 *
 * @private
 * @param {Geometry|Feature<any>} geojson Base projection on center of GeoJSON
 * @returns {GeoProjection} D3 Geo Transverse Mercator Projection
 */
function defineProjection(geojson) {
    var coords = center(geojson).geometry.coordinates.reverse();
    var rotate = coords.map(function (coord) { return -coord; });
    return d3Geo.geoTransverseMercator()
        .center(coords)
        .rotate(rotate)
        .scale(helpers.earthRadius);
}

module.exports = buffer;
module.exports.default = buffer;

},{"@turf/bbox":3,"@turf/center":8,"@turf/helpers":7,"@turf/meta":34,"@turf/projection":40,"d3-geo":50,"turf-jsts":67}],7:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

/**
 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
 */
var earthRadius = 6371008.8;

/**
 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
 */
var factors = {
    meters: earthRadius,
    metres: earthRadius,
    millimeters: earthRadius * 1000,
    millimetres: earthRadius * 1000,
    centimeters: earthRadius * 100,
    centimetres: earthRadius * 100,
    kilometers: earthRadius / 1000,
    kilometres: earthRadius / 1000,
    miles: earthRadius / 1609.344,
    nauticalmiles: earthRadius / 1852,
    inches: earthRadius * 39.370,
    yards: earthRadius / 1.0936,
    feet: earthRadius * 3.28084,
    radians: 1,
    degrees: earthRadius / 111325,
};

/**
 * Units of measurement factors based on 1 meter.
 */
var unitsFactors = {
    meters: 1,
    metres: 1,
    millimeters: 1000,
    millimetres: 1000,
    centimeters: 100,
    centimetres: 100,
    kilometers: 1 / 1000,
    kilometres: 1 / 1000,
    miles: 1 / 1609.344,
    nauticalmiles: 1 / 1852,
    inches: 39.370,
    yards: 1 / 1.0936,
    feet: 3.28084,
    radians: 1 / earthRadius,
    degrees: 1 / 111325,
};

/**
 * Area of measurement factors based on 1 square meter.
 */
var areaFactors = {
    meters: 1,
    metres: 1,
    millimeters: 1000000,
    millimetres: 1000000,
    centimeters: 10000,
    centimetres: 10000,
    kilometers: 0.000001,
    kilometres: 0.000001,
    acres: 0.000247105,
    miles: 3.86e-7,
    yards: 1.195990046,
    feet: 10.763910417,
    inches: 1550.003100006
};

/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature} a GeoJSON Feature
 * @example
 * var geometry = {
 *   "type": "Point",
 *   "coordinates": [110, 50]
 * };
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geometry, properties, options) {
    // Optional Parameters
    options = options || {};
    if (!isObject(options)) throw new Error('options is invalid');
    var bbox = options.bbox;
    var id = options.id;

    // Validation
    if (geometry === undefined) throw new Error('geometry is required');
    if (properties && properties.constructor !== Object) throw new Error('properties must be an Object');
    if (bbox) validateBBox(bbox);
    if (id) validateId(id);

    // Main
    var feat = {type: 'Feature'};
    if (id) feat.id = id;
    if (bbox) feat.bbox = bbox;
    feat.properties = properties || {};
    feat.geometry = geometry;
    return feat;
}

/**
 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
 * For GeometryCollection type use `helpers.geometryCollection`
 *
 * @name geometry
 * @param {string} type Geometry Type
 * @param {Array<number>} coordinates Coordinates
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Geometry
 * @returns {Geometry} a GeoJSON Geometry
 * @example
 * var type = 'Point';
 * var coordinates = [110, 50];
 *
 * var geometry = turf.geometry(type, coordinates);
 *
 * //=geometry
 */
function geometry(type, coordinates, options) {
    // Optional Parameters
    options = options || {};
    if (!isObject(options)) throw new Error('options is invalid');
    var bbox = options.bbox;

    // Validation
    if (!type) throw new Error('type is required');
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');
    if (bbox) validateBBox(bbox);

    // Main
    var geom;
    switch (type) {
    case 'Point': geom = point(coordinates).geometry; break;
    case 'LineString': geom = lineString(coordinates).geometry; break;
    case 'Polygon': geom = polygon(coordinates).geometry; break;
    case 'MultiPoint': geom = multiPoint(coordinates).geometry; break;
    case 'MultiLineString': geom = multiLineString(coordinates).geometry; break;
    case 'MultiPolygon': geom = multiPolygon(coordinates).geometry; break;
    default: throw new Error(type + ' is invalid');
    }
    if (bbox) geom.bbox = bbox;
    return geom;
}

/**
 * Creates a {@link Point} {@link Feature} from a Position.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Point>} a Point feature
 * @example
 * var point = turf.point([-75.343, 39.984]);
 *
 * //=point
 */
function point(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');
    if (coordinates.length < 2) throw new Error('coordinates must be at least 2 numbers long');
    if (!isNumber(coordinates[0]) || !isNumber(coordinates[1])) throw new Error('coordinates must contain numbers');

    return feature({
        type: 'Point',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
 *
 * @name points
 * @param {Array<Array<number>>} coordinates an array of Points
 * @param {Object} [properties={}] Translate these properties to each Feature
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Point>} Point Feature
 * @example
 * var points = turf.points([
 *   [-75, 39],
 *   [-80, 45],
 *   [-78, 50]
 * ]);
 *
 * //=points
 */
function points(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');

    return featureCollection(coordinates.map(function (coords) {
        return point(coords, properties);
    }), options);
}

/**
 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Polygon>} Polygon Feature
 * @example
 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
 *
 * //=polygon
 */
function polygon(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');

    for (var i = 0; i < coordinates.length; i++) {
        var ring = coordinates[i];
        if (ring.length < 4) {
            throw new Error('Each LinearRing of a Polygon must have 4 or more Positions.');
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            // Check if first point of Polygon contains two numbers
            if (i === 0 && j === 0 && !isNumber(ring[0][0]) || !isNumber(ring[0][1])) throw new Error('coordinates must contain numbers');
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error('First and last Position are not equivalent.');
            }
        }
    }

    return feature({
        type: 'Polygon',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
 *
 * @name polygons
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
 * @example
 * var polygons = turf.polygons([
 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
 * ]);
 *
 * //=polygons
 */
function polygons(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');

    return featureCollection(coordinates.map(function (coords) {
        return polygon(coords, properties);
    }), options);
}

/**
 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<LineString>} LineString Feature
 * @example
 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
 *
 * //=linestring1
 * //=linestring2
 */
function lineString(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');
    if (coordinates.length < 2) throw new Error('coordinates must be an array of two or more positions');
    // Check if first point of LineString contains two numbers
    if (!isNumber(coordinates[0][1]) || !isNumber(coordinates[0][1])) throw new Error('coordinates must contain numbers');

    return feature({
        type: 'LineString',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
 *
 * @name lineStrings
 * @param {Array<Array<number>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
 * @example
 * var linestrings = turf.lineStrings([
 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
 * ]);
 *
 * //=linestrings
 */
function lineStrings(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');
    if (!Array.isArray(coordinates)) throw new Error('coordinates must be an Array');

    return featureCollection(coordinates.map(function (coords) {
        return lineString(coords, properties);
    }), options);
}

/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {FeatureCollection} FeatureCollection of Features
 * @example
 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
 *
 * var collection = turf.featureCollection([
 *   locationA,
 *   locationB,
 *   locationC
 * ]);
 *
 * //=collection
 */
function featureCollection(features, options) {
    // Optional Parameters
    options = options || {};
    if (!isObject(options)) throw new Error('options is invalid');
    var bbox = options.bbox;
    var id = options.id;

    // Validation
    if (!features) throw new Error('No features passed');
    if (!Array.isArray(features)) throw new Error('features must be an Array');
    if (bbox) validateBBox(bbox);
    if (id) validateId(id);

    // Main
    var fc = {type: 'FeatureCollection'};
    if (id) fc.id = id;
    if (bbox) fc.bbox = bbox;
    fc.features = features;
    return fc;
}

/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 */
function multiLineString(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');

    return feature({
        type: 'MultiLineString',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 */
function multiPoint(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');

    return feature({
        type: 'MultiPoint',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
function multiPolygon(coordinates, properties, options) {
    if (!coordinates) throw new Error('coordinates is required');

    return feature({
        type: 'MultiPolygon',
        coordinates: coordinates
    }, properties, options);
}

/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = {
 *     "type": "Point",
 *       "coordinates": [100, 0]
 *     };
 * var line = {
 *     "type": "LineString",
 *     "coordinates": [ [101, 0], [102, 1] ]
 *   };
 * var collection = turf.geometryCollection([pt, line]);
 *
 * //=collection
 */
function geometryCollection(geometries, properties, options) {
    if (!geometries) throw new Error('geometries is required');
    if (!Array.isArray(geometries)) throw new Error('geometries must be an Array');

    return feature({
        type: 'GeometryCollection',
        geometries: geometries
    }, properties, options);
}

/**
 * Round number to precision
 *
 * @param {number} num Number
 * @param {number} [precision=0] Precision
 * @returns {number} rounded number
 * @example
 * turf.round(120.4321)
 * //=120
 *
 * turf.round(120.4321, 2)
 * //=120.43
 */
function round(num, precision) {
    if (num === undefined || num === null || isNaN(num)) throw new Error('num is required');
    if (precision && !(precision >= 0)) throw new Error('precision must be a positive number');
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(num * multiplier) / multiplier;
}

/**
 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name radiansToLength
 * @param {number} radians in radians across the sphere
 * @param {string} [units='kilometers'] can be degrees, radians, miles, or kilometers inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} distance
 */
function radiansToLength(radians, units) {
    if (radians === undefined || radians === null) throw new Error('radians is required');

    if (units && typeof units !== 'string') throw new Error('units must be a string');
    var factor = factors[units || 'kilometers'];
    if (!factor) throw new Error(units + ' units is invalid');
    return radians * factor;
}

/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name lengthToRadians
 * @param {number} distance in real units
 * @param {string} [units='kilometers'] can be degrees, radians, miles, or kilometers inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} radians
 */
function lengthToRadians(distance, units) {
    if (distance === undefined || distance === null) throw new Error('distance is required');

    if (units && typeof units !== 'string') throw new Error('units must be a string');
    var factor = factors[units || 'kilometers'];
    if (!factor) throw new Error(units + ' units is invalid');
    return distance / factor;
}

/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
 *
 * @name lengthToDegrees
 * @param {number} distance in real units
 * @param {string} [units='kilometers'] can be degrees, radians, miles, or kilometers inches, yards, metres, meters, kilometres, kilometers.
 * @returns {number} degrees
 */
function lengthToDegrees(distance, units) {
    return radiansToDegrees(lengthToRadians(distance, units));
}

/**
 * Converts any bearing angle from the north line direction (positive clockwise)
 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
 *
 * @name bearingToAzimuth
 * @param {number} bearing angle, between -180 and +180 degrees
 * @returns {number} angle between 0 and 360 degrees
 */
function bearingToAzimuth(bearing) {
    if (bearing === null || bearing === undefined) throw new Error('bearing is required');

    var angle = bearing % 360;
    if (angle < 0) angle += 360;
    return angle;
}

/**
 * Converts an angle in radians to degrees
 *
 * @name radiansToDegrees
 * @param {number} radians angle in radians
 * @returns {number} degrees between 0 and 360 degrees
 */
function radiansToDegrees(radians) {
    if (radians === null || radians === undefined) throw new Error('radians is required');

    var degrees = radians % (2 * Math.PI);
    return degrees * 180 / Math.PI;
}

/**
 * Converts an angle in degrees to radians
 *
 * @name degreesToRadians
 * @param {number} degrees angle between 0 and 360 degrees
 * @returns {number} angle in radians
 */
function degreesToRadians(degrees) {
    if (degrees === null || degrees === undefined) throw new Error('degrees is required');

    var radians = degrees % 360;
    return radians * Math.PI / 180;
}

/**
 * Converts a length to the requested unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @param {number} length to be converted
 * @param {string} originalUnit of the length
 * @param {string} [finalUnit='kilometers'] returned unit
 * @returns {number} the converted length
 */
function convertLength(length, originalUnit, finalUnit) {
    if (length === null || length === undefined) throw new Error('length is required');
    if (!(length >= 0)) throw new Error('length must be a positive number');

    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit || 'kilometers');
}

/**
 * Converts a area to the requested unit.
 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches
 * @param {number} area to be converted
 * @param {string} [originalUnit='meters'] of the distance
 * @param {string} [finalUnit='kilometers'] returned unit
 * @returns {number} the converted distance
 */
function convertArea(area, originalUnit, finalUnit) {
    if (area === null || area === undefined) throw new Error('area is required');
    if (!(area >= 0)) throw new Error('area must be a positive number');

    var startFactor = areaFactors[originalUnit || 'meters'];
    if (!startFactor) throw new Error('invalid original units');

    var finalFactor = areaFactors[finalUnit || 'kilometers'];
    if (!finalFactor) throw new Error('invalid final units');

    return (area / startFactor) * finalFactor;
}

/**
 * isNumber
 *
 * @param {*} num Number to validate
 * @returns {boolean} true/false
 * @example
 * turf.isNumber(123)
 * //=true
 * turf.isNumber('foo')
 * //=false
 */
function isNumber(num) {
    return !isNaN(num) && num !== null && !Array.isArray(num);
}

/**
 * isObject
 *
 * @param {*} input variable to validate
 * @returns {boolean} true/false
 * @example
 * turf.isObject({elevation: 10})
 * //=true
 * turf.isObject('foo')
 * //=false
 */
function isObject(input) {
    return (!!input) && (input.constructor === Object);
}

/**
 * Validate BBox
 *
 * @private
 * @param {Array<number>} bbox BBox to validate
 * @returns {void}
 * @throws Error if BBox is not valid
 * @example
 * validateBBox([-180, -40, 110, 50])
 * //=OK
 * validateBBox([-180, -40])
 * //=Error
 * validateBBox('Foo')
 * //=Error
 * validateBBox(5)
 * //=Error
 * validateBBox(null)
 * //=Error
 * validateBBox(undefined)
 * //=Error
 */
function validateBBox(bbox) {
    if (!bbox) throw new Error('bbox is required');
    if (!Array.isArray(bbox)) throw new Error('bbox must be an Array');
    if (bbox.length !== 4 && bbox.length !== 6) throw new Error('bbox must be an Array of 4 or 6 numbers');
    bbox.forEach(function (num) {
        if (!isNumber(num)) throw new Error('bbox must only contain numbers');
    });
}

/**
 * Validate Id
 *
 * @private
 * @param {string|number} id Id to validate
 * @returns {void}
 * @throws Error if Id is not valid
 * @example
 * validateId([-180, -40, 110, 50])
 * //=Error
 * validateId([-180, -40])
 * //=Error
 * validateId('Foo')
 * //=OK
 * validateId(5)
 * //=OK
 * validateId(null)
 * //=Error
 * validateId(undefined)
 * //=Error
 */
function validateId(id) {
    if (!id) throw new Error('id is required');
    if (['string', 'number'].indexOf(typeof id) === -1) throw new Error('id must be a number or a string');
}

// Deprecated methods
function radians2degrees() {
    throw new Error('method has been renamed to `radiansToDegrees`');
}

function degrees2radians() {
    throw new Error('method has been renamed to `degreesToRadians`');
}

function distanceToDegrees() {
    throw new Error('method has been renamed to `lengthToDegrees`');
}

function distanceToRadians() {
    throw new Error('method has been renamed to `lengthToRadians`');
}

function radiansToDistance() {
    throw new Error('method has been renamed to `radiansToLength`');
}

function bearingToAngle() {
    throw new Error('method has been renamed to `bearingToAzimuth`');
}

function convertDistance() {
    throw new Error('method has been renamed to `convertLength`');
}

exports.earthRadius = earthRadius;
exports.factors = factors;
exports.unitsFactors = unitsFactors;
exports.areaFactors = areaFactors;
exports.feature = feature;
exports.geometry = geometry;
exports.point = point;
exports.points = points;
exports.polygon = polygon;
exports.polygons = polygons;
exports.lineString = lineString;
exports.lineStrings = lineStrings;
exports.featureCollection = featureCollection;
exports.multiLineString = multiLineString;
exports.multiPoint = multiPoint;
exports.multiPolygon = multiPolygon;
exports.geometryCollection = geometryCollection;
exports.round = round;
exports.radiansToLength = radiansToLength;
exports.lengthToRadians = lengthToRadians;
exports.lengthToDegrees = lengthToDegrees;
exports.bearingToAzimuth = bearingToAzimuth;
exports.radiansToDegrees = radiansToDegrees;
exports.degreesToRadians = degreesToRadians;
exports.convertLength = convertLength;
exports.convertArea = convertArea;
exports.isNumber = isNumber;
exports.isObject = isObject;
exports.validateBBox = validateBBox;
exports.validateId = validateId;
exports.radians2degrees = radians2degrees;
exports.degrees2radians = degrees2radians;
exports.distanceToDegrees = distanceToDegrees;
exports.distanceToRadians = distanceToRadians;
exports.radiansToDistance = radiansToDistance;
exports.bearingToAngle = bearingToAngle;
exports.convertDistance = convertDistance;

},{}],8:[function(require,module,exports){
'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var bbox = _interopDefault(require('@turf/bbox'));
var helpers = require('@turf/helpers');

/**
 * Takes a {@link Feature} or {@link FeatureCollection} and returns the absolute center point of all features.
 *
 * @name center
 * @param {GeoJSON} geojson GeoJSON to be centered
 * @param {Object} [options={}] Optional parameters
 * @param {Object} [options.properties={}] an Object that is used as the {@link Feature}'s properties
 * @returns {Feature<Point>} a Point feature at the absolute center point of all input features
 * @example
 * var features = turf.featureCollection([
 *   turf.point( [-97.522259, 35.4691]),
 *   turf.point( [-97.502754, 35.463455]),
 *   turf.point( [-97.508269, 35.463245])
 * ]);
 *
 * var center = turf.center(features);
 *
 * //addToMap
 * var addToMap = [features, center]
 * center.properties['marker-size'] = 'large';
 * center.properties['marker-color'] = '#000';
 */
function center(geojson, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var properties = options.properties;

    // Input validation
    if (!geojson) throw new Error('geojson is required');

    var ext = bbox(geojson);
    var x = (ext[0] + ext[2]) / 2;
    var y = (ext[1] + ext[3]) / 2;
    return helpers.point([x, y], properties);
}

module.exports = center;
module.exports.default = center;

},{"@turf/bbox":3,"@turf/helpers":9}],9:[function(require,module,exports){
arguments[4][7][0].apply(exports,arguments)
},{"dup":7}],10:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var meta_1 = require("@turf/meta");
var helpers_1 = require("@turf/helpers");
/**
 * Takes one or more features and calculates the centroid using the mean of all vertices.
 * This lessens the effect of small islands and artifacts when calculating the centroid of a set of polygons.
 *
 * @name centroid
 * @param {GeoJSON} geojson GeoJSON to be centered
 * @param {Object} [options={}] Optional Parameters
 * @param {Object} [options.properties={}] an Object that is used as the {@link Feature}'s properties
 * @returns {Feature<Point>} the centroid of the input features
 * @example
 * var polygon = turf.polygon([[[-81, 41], [-88, 36], [-84, 31], [-80, 33], [-77, 39], [-81, 41]]]);
 *
 * var centroid = turf.centroid(polygon);
 *
 * //addToMap
 * var addToMap = [polygon, centroid]
 */
function centroid(geojson, options) {
    if (options === void 0) { options = {}; }
    var xSum = 0;
    var ySum = 0;
    var len = 0;
    meta_1.coordEach(geojson, function (coord) {
        xSum += coord[0];
        ySum += coord[1];
        len++;
    });
    return helpers_1.point([xSum / len, ySum / len], options.properties);
}
exports.default = centroid;

},{"@turf/helpers":15,"@turf/meta":11}],11:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

var helpers = require('@turf/helpers');

/**
 * Callback for coordEach
 *
 * @callback coordEachCallback
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
 *
 * @name coordEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function coordEach(geojson, callback, excludeWrapCoord) {
    // Handles null Geometry -- Skips this GeoJSON
    if (geojson === null) return;
    var j, k, l, geometry, stopG, coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        coordIndex = 0,
        isGeometryCollection,
        type = geojson.type,
        isFeatureCollection = type === 'FeatureCollection',
        isFeature = type === 'Feature',
        stop = isFeatureCollection ? geojson.features.length : 1;

    // This logic may look a little weird. The reason why it is that way
    // is because it's trying to be fast. GeoJSON supports multiple kinds
    // of objects at its root: FeatureCollection, Features, Geometries.
    // This function has the responsibility of handling all of them, and that
    // means that some of the `for` loops you see below actually just don't apply
    // to certain inputs. For instance, if you give this just a
    // Point geometry, then both loops are short-circuited and all we do
    // is gradually rename the input until it's called 'geometry'.
    //
    // This also aims to allocate as few resources as possible: just a
    // few numbers and booleans, rather than any temporary arrays as would
    // be required with the normalization approach.
    for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = (isFeatureCollection ? geojson.features[featureIndex].geometry :
            (isFeature ? geojson.geometry : geojson));
        isGeometryCollection = (geometryMaybeCollection) ? geometryMaybeCollection.type === 'GeometryCollection' : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
            var multiFeatureIndex = 0;
            var geometryIndex = 0;
            geometry = isGeometryCollection ?
                geometryMaybeCollection.geometries[geomIndex] : geometryMaybeCollection;

            // Handles null Geometry -- Skips this geometry
            if (geometry === null) continue;
            coords = geometry.coordinates;
            var geomType = geometry.type;

            wrapShrink = (excludeWrapCoord && (geomType === 'Polygon' || geomType === 'MultiPolygon')) ? 1 : 0;

            switch (geomType) {
            case null:
                break;
            case 'Point':
                if (callback(coords, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                coordIndex++;
                multiFeatureIndex++;
                break;
            case 'LineString':
            case 'MultiPoint':
                for (j = 0; j < coords.length; j++) {
                    if (callback(coords[j], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                    coordIndex++;
                    if (geomType === 'MultiPoint') multiFeatureIndex++;
                }
                if (geomType === 'LineString') multiFeatureIndex++;
                break;
            case 'Polygon':
            case 'MultiLineString':
                for (j = 0; j < coords.length; j++) {
                    for (k = 0; k < coords[j].length - wrapShrink; k++) {
                        if (callback(coords[j][k], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                        coordIndex++;
                    }
                    if (geomType === 'MultiLineString') multiFeatureIndex++;
                    if (geomType === 'Polygon') geometryIndex++;
                }
                if (geomType === 'Polygon') multiFeatureIndex++;
                break;
            case 'MultiPolygon':
                for (j = 0; j < coords.length; j++) {
                    geometryIndex = 0;
                    for (k = 0; k < coords[j].length; k++) {
                        for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                            if (callback(coords[j][k][l], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                            coordIndex++;
                        }
                        geometryIndex++;
                    }
                    multiFeatureIndex++;
                }
                break;
            case 'GeometryCollection':
                for (j = 0; j < geometry.geometries.length; j++)
                    if (coordEach(geometry.geometries[j], callback, excludeWrapCoord) === false) return false;
                break;
            default:
                throw new Error('Unknown Geometry Type');
            }
        }
    }
}

/**
 * Callback for coordReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback coordReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Reduce coordinates in any GeoJSON object, similar to Array.reduce()
 *
 * @name coordReduce
 * @param {FeatureCollection|Geometry|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentCoord, coordIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordReduce(features, function (previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentCoord;
 * });
 */
function coordReduce(geojson, callback, initialValue, excludeWrapCoord) {
    var previousValue = initialValue;
    coordEach(geojson, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
        if (coordIndex === 0 && initialValue === undefined) previousValue = currentCoord;
        else previousValue = callback(previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex);
    }, excludeWrapCoord);
    return previousValue;
}

/**
 * Callback for propEach
 *
 * @callback propEachCallback
 * @param {Object} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over properties in any GeoJSON object, similar to Array.forEach()
 *
 * @name propEach
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentProperties, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propEach(features, function (currentProperties, featureIndex) {
 *   //=currentProperties
 *   //=featureIndex
 * });
 */
function propEach(geojson, callback) {
    var i;
    switch (geojson.type) {
    case 'FeatureCollection':
        for (i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i].properties, i) === false) break;
        }
        break;
    case 'Feature':
        callback(geojson.properties, 0);
        break;
    }
}


/**
 * Callback for propReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback propReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {*} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce properties in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all properties is unnecessary.
 *
 * @name propReduce
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentProperties, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propReduce(features, function (previousValue, currentProperties, featureIndex) {
 *   //=previousValue
 *   //=currentProperties
 *   //=featureIndex
 *   return currentProperties
 * });
 */
function propReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    propEach(geojson, function (currentProperties, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentProperties;
        else previousValue = callback(previousValue, currentProperties, featureIndex);
    });
    return previousValue;
}

/**
 * Callback for featureEach
 *
 * @callback featureEachCallback
 * @param {Feature<any>} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name featureEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.featureEach(features, function (currentFeature, featureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 * });
 */
function featureEach(geojson, callback) {
    if (geojson.type === 'Feature') {
        callback(geojson, 0);
    } else if (geojson.type === 'FeatureCollection') {
        for (var i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i], i) === false) break;
        }
    }
}

/**
 * Callback for featureReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback featureReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name featureReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.featureReduce(features, function (previousValue, currentFeature, featureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   return currentFeature
 * });
 */
function featureReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    featureEach(geojson, function (currentFeature, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentFeature;
        else previousValue = callback(previousValue, currentFeature, featureIndex);
    });
    return previousValue;
}

/**
 * Get all coordinates from any GeoJSON object.
 *
 * @name coordAll
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @returns {Array<Array<number>>} coordinate position array
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * var coords = turf.coordAll(features);
 * //= [[26, 37], [36, 53]]
 */
function coordAll(geojson) {
    var coords = [];
    coordEach(geojson, function (coord) {
        coords.push(coord);
    });
    return coords;
}

/**
 * Callback for geomEach
 *
 * @callback geomEachCallback
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
 *
 * @name geomEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomEach(features, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 * });
 */
function geomEach(geojson, callback) {
    var i, j, g, geometry, stopG,
        geometryMaybeCollection,
        isGeometryCollection,
        featureProperties,
        featureBBox,
        featureId,
        featureIndex = 0,
        isFeatureCollection = geojson.type === 'FeatureCollection',
        isFeature = geojson.type === 'Feature',
        stop = isFeatureCollection ? geojson.features.length : 1;

    // This logic may look a little weird. The reason why it is that way
    // is because it's trying to be fast. GeoJSON supports multiple kinds
    // of objects at its root: FeatureCollection, Features, Geometries.
    // This function has the responsibility of handling all of them, and that
    // means that some of the `for` loops you see below actually just don't apply
    // to certain inputs. For instance, if you give this just a
    // Point geometry, then both loops are short-circuited and all we do
    // is gradually rename the input until it's called 'geometry'.
    //
    // This also aims to allocate as few resources as possible: just a
    // few numbers and booleans, rather than any temporary arrays as would
    // be required with the normalization approach.
    for (i = 0; i < stop; i++) {

        geometryMaybeCollection = (isFeatureCollection ? geojson.features[i].geometry :
            (isFeature ? geojson.geometry : geojson));
        featureProperties = (isFeatureCollection ? geojson.features[i].properties :
            (isFeature ? geojson.properties : {}));
        featureBBox = (isFeatureCollection ? geojson.features[i].bbox :
            (isFeature ? geojson.bbox : undefined));
        featureId = (isFeatureCollection ? geojson.features[i].id :
            (isFeature ? geojson.id : undefined));
        isGeometryCollection = (geometryMaybeCollection) ? geometryMaybeCollection.type === 'GeometryCollection' : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (g = 0; g < stopG; g++) {
            geometry = isGeometryCollection ?
                geometryMaybeCollection.geometries[g] : geometryMaybeCollection;

            // Handle null Geometry
            if (geometry === null) {
                if (callback(null, featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                continue;
            }
            switch (geometry.type) {
            case 'Point':
            case 'LineString':
            case 'MultiPoint':
            case 'Polygon':
            case 'MultiLineString':
            case 'MultiPolygon': {
                if (callback(geometry, featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                break;
            }
            case 'GeometryCollection': {
                for (j = 0; j < geometry.geometries.length; j++) {
                    if (callback(geometry.geometries[j], featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                }
                break;
            }
            default:
                throw new Error('Unknown Geometry Type');
            }
        }
        // Only increase `featureIndex` per each feature
        featureIndex++;
    }
}

/**
 * Callback for geomReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback geomReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Reduce geometry in any GeoJSON object, similar to Array.reduce().
 *
 * @name geomReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomReduce(features, function (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=previousValue
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 *   return currentGeometry
 * });
 */
function geomReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    geomEach(geojson, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentGeometry;
        else previousValue = callback(previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId);
    });
    return previousValue;
}

/**
 * Callback for flattenEach
 *
 * @callback flattenEachCallback
 * @param {Feature} currentFeature The current flattened feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Iterate over flattened features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name flattenEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex, multiFeatureIndex)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenEach(features, function (currentFeature, featureIndex, multiFeatureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 * });
 */
function flattenEach(geojson, callback) {
    geomEach(geojson, function (geometry, featureIndex, properties, bbox, id) {
        // Callback for single geometry
        var type = (geometry === null) ? null : geometry.type;
        switch (type) {
        case null:
        case 'Point':
        case 'LineString':
        case 'Polygon':
            if (callback(helpers.feature(geometry, properties, {bbox: bbox, id: id}), featureIndex, 0) === false) return false;
            return;
        }

        var geomType;

        // Callback for multi-geometry
        switch (type) {
        case 'MultiPoint':
            geomType = 'Point';
            break;
        case 'MultiLineString':
            geomType = 'LineString';
            break;
        case 'MultiPolygon':
            geomType = 'Polygon';
            break;
        }

        for (var multiFeatureIndex = 0; multiFeatureIndex < geometry.coordinates.length; multiFeatureIndex++) {
            var coordinate = geometry.coordinates[multiFeatureIndex];
            var geom = {
                type: geomType,
                coordinates: coordinate
            };
            if (callback(helpers.feature(geom, properties), featureIndex, multiFeatureIndex) === false) return false;
        }
    });
}

/**
 * Callback for flattenReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback flattenReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Reduce flattened features in any GeoJSON object, similar to Array.reduce().
 *
 * @name flattenReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex, multiFeatureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenReduce(features, function (previousValue, currentFeature, featureIndex, multiFeatureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   return currentFeature
 * });
 */
function flattenReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    flattenEach(geojson, function (currentFeature, featureIndex, multiFeatureIndex) {
        if (featureIndex === 0 && multiFeatureIndex === 0 && initialValue === undefined) previousValue = currentFeature;
        else previousValue = callback(previousValue, currentFeature, featureIndex, multiFeatureIndex);
    });
    return previousValue;
}

/**
 * Callback for segmentEach
 *
 * @callback segmentEachCallback
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 * @returns {void}
 */

/**
 * Iterate over 2-vertex line segment in any GeoJSON object, similar to Array.forEach()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex)
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentEach(polygon, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //=currentSegment
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   //=segmentIndex
 * });
 *
 * // Calculate the total number of segments
 * var total = 0;
 * turf.segmentEach(polygon, function () {
 *     total++;
 * });
 */
function segmentEach(geojson, callback) {
    flattenEach(geojson, function (feature, featureIndex, multiFeatureIndex) {
        var segmentIndex = 0;

        // Exclude null Geometries
        if (!feature.geometry) return;
        // (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
        var type = feature.geometry.type;
        if (type === 'Point' || type === 'MultiPoint') return;

        // Generate 2-vertex line segments
        var previousCoords;
        var previousFeatureIndex = 0;
        var previousMultiIndex = 0;
        var prevGeomIndex = 0;
        if (coordEach(feature, function (currentCoord, coordIndex, featureIndexCoord, multiPartIndexCoord, geometryIndex) {
            // Simulating a meta.coordReduce() since `reduce` operations cannot be stopped by returning `false`
            if (previousCoords === undefined || featureIndex > previousFeatureIndex || multiPartIndexCoord > previousMultiIndex || geometryIndex > prevGeomIndex) {
                previousCoords = currentCoord;
                previousFeatureIndex = featureIndex;
                previousMultiIndex = multiPartIndexCoord;
                prevGeomIndex = geometryIndex;
                segmentIndex = 0;
                return;
            }
            var currentSegment = helpers.lineString([previousCoords, currentCoord], feature.properties);
            if (callback(currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) === false) return false;
            segmentIndex++;
            previousCoords = currentCoord;
        }) === false) return false;
    });
}

/**
 * Callback for segmentReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback segmentReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 */

/**
 * Reduce 2-vertex line segment in any GeoJSON object, similar to Array.reduce()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (previousValue, currentSegment, currentIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentReduce(polygon, function (previousSegment, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //= previousSegment
 *   //= currentSegment
 *   //= featureIndex
 *   //= multiFeatureIndex
 *   //= geometryIndex
 *   //= segmentInex
 *   return currentSegment
 * });
 *
 * // Calculate the total number of segments
 * var initialValue = 0
 * var total = turf.segmentReduce(polygon, function (previousValue) {
 *     previousValue++;
 *     return previousValue;
 * }, initialValue);
 */
function segmentReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    var started = false;
    segmentEach(geojson, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
        if (started === false && initialValue === undefined) previousValue = currentSegment;
        else previousValue = callback(previousValue, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex);
        started = true;
    });
    return previousValue;
}

/**
 * Callback for lineEach
 *
 * @callback lineEachCallback
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Iterate over line or ring coordinates in LineString, Polygon, MultiLineString, MultiPolygon Features or Geometries,
 * similar to Array.forEach.
 *
 * @name lineEach
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @example
 * var multiLine = turf.multiLineString([
 *   [[26, 37], [35, 45]],
 *   [[36, 53], [38, 50], [41, 55]]
 * ]);
 *
 * turf.lineEach(multiLine, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function lineEach(geojson, callback) {
    // validation
    if (!geojson) throw new Error('geojson is required');

    flattenEach(geojson, function (feature, featureIndex, multiFeatureIndex) {
        if (feature.geometry === null) return;
        var type = feature.geometry.type;
        var coords = feature.geometry.coordinates;
        switch (type) {
        case 'LineString':
            if (callback(feature, featureIndex, multiFeatureIndex, 0, 0) === false) return false;
            break;
        case 'Polygon':
            for (var geometryIndex = 0; geometryIndex < coords.length; geometryIndex++) {
                if (callback(helpers.lineString(coords[geometryIndex], feature.properties), featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
            }
            break;
        }
    });
}

/**
 * Callback for lineReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback lineReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name lineReduce
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var multiPoly = turf.multiPolygon([
 *   turf.polygon([[[12,48],[2,41],[24,38],[12,48]], [[9,44],[13,41],[13,45],[9,44]]]),
 *   turf.polygon([[[5, 5], [0, 0], [2, 2], [4, 4], [5, 5]]])
 * ]);
 *
 * turf.lineReduce(multiPoly, function (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentLine
 * });
 */
function lineReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    lineEach(geojson, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentLine;
        else previousValue = callback(previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex);
    });
    return previousValue;
}

/**
 * Finds a particular 2-vertex LineString Segment from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 * Point & MultiPoint will always return null.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.segmentIndex=0] Segment Index
 * @param {Object} [options.properties={}] Translate Properties to output LineString
 * @param {BBox} [options.bbox={}] Translate BBox to output LineString
 * @param {number|string} [options.id={}] Translate Id to output LineString
 * @returns {Feature<LineString>} 2-vertex GeoJSON Feature LineString
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findSegment(multiLine);
 * // => Feature<LineString<[[10, 10], [50, 30]]>>
 *
 * // First Segment of 2nd Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: 1});
 * // => Feature<LineString<[[-10, -10], [-50, -30]]>>
 *
 * // Last Segment of Last Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: -1, segmentIndex: -1});
 * // => Feature<LineString<[[-50, -30], [-30, -40]]>>
 */
function findSegment(geojson, options) {
    // Optional Parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var featureIndex = options.featureIndex || 0;
    var multiFeatureIndex = options.multiFeatureIndex || 0;
    var geometryIndex = options.geometryIndex || 0;
    var segmentIndex = options.segmentIndex || 0;

    // Find FeatureIndex
    var properties = options.properties;
    var geometry;

    switch (geojson.type) {
    case 'FeatureCollection':
        if (featureIndex < 0) featureIndex = geojson.features.length + featureIndex;
        properties = properties || geojson.features[featureIndex].properties;
        geometry = geojson.features[featureIndex].geometry;
        break;
    case 'Feature':
        properties = properties || geojson.properties;
        geometry = geojson.geometry;
        break;
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
    case 'Polygon':
    case 'MultiLineString':
    case 'MultiPolygon':
        geometry = geojson;
        break;
    default:
        throw new Error('geojson is invalid');
    }

    // Find SegmentIndex
    if (geometry === null) return null;
    var coords = geometry.coordinates;
    switch (geometry.type) {
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
        if (segmentIndex < 0) segmentIndex = coords.length + segmentIndex - 1;
        return helpers.lineString([coords[segmentIndex], coords[segmentIndex + 1]], properties, options);
    case 'Polygon':
        if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
        if (segmentIndex < 0) segmentIndex = coords[geometryIndex].length + segmentIndex - 1;
        return helpers.lineString([coords[geometryIndex][segmentIndex], coords[geometryIndex][segmentIndex + 1]], properties, options);
    case 'MultiLineString':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (segmentIndex < 0) segmentIndex = coords[multiFeatureIndex].length + segmentIndex - 1;
        return helpers.lineString([coords[multiFeatureIndex][segmentIndex], coords[multiFeatureIndex][segmentIndex + 1]], properties, options);
    case 'MultiPolygon':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (geometryIndex < 0) geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
        if (segmentIndex < 0) segmentIndex = coords[multiFeatureIndex][geometryIndex].length - segmentIndex - 1;
        return helpers.lineString([coords[multiFeatureIndex][geometryIndex][segmentIndex], coords[multiFeatureIndex][geometryIndex][segmentIndex + 1]], properties, options);
    }
    throw new Error('geojson is invalid');
}

/**
 * Finds a particular Point from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.coordIndex=0] Coord Index
 * @param {Object} [options.properties={}] Translate Properties to output Point
 * @param {BBox} [options.bbox={}] Translate BBox to output Point
 * @param {number|string} [options.id={}] Translate Id to output Point
 * @returns {Feature<Point>} 2-vertex GeoJSON Feature Point
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findPoint(multiLine);
 * // => Feature<Point<[10, 10]>>
 *
 * // First Segment of the 2nd Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: 1});
 * // => Feature<Point<[-10, -10]>>
 *
 * // Last Segment of last Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: -1, coordIndex: -1});
 * // => Feature<Point<[-30, -40]>>
 */
function findPoint(geojson, options) {
    // Optional Parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var featureIndex = options.featureIndex || 0;
    var multiFeatureIndex = options.multiFeatureIndex || 0;
    var geometryIndex = options.geometryIndex || 0;
    var coordIndex = options.coordIndex || 0;

    // Find FeatureIndex
    var properties = options.properties;
    var geometry;

    switch (geojson.type) {
    case 'FeatureCollection':
        if (featureIndex < 0) featureIndex = geojson.features.length + featureIndex;
        properties = properties || geojson.features[featureIndex].properties;
        geometry = geojson.features[featureIndex].geometry;
        break;
    case 'Feature':
        properties = properties || geojson.properties;
        geometry = geojson.geometry;
        break;
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
    case 'Polygon':
    case 'MultiLineString':
    case 'MultiPolygon':
        geometry = geojson;
        break;
    default:
        throw new Error('geojson is invalid');
    }

    // Find Coord Index
    if (geometry === null) return null;
    var coords = geometry.coordinates;
    switch (geometry.type) {
    case 'Point':
        return helpers.point(coords, properties, options);
    case 'MultiPoint':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        return helpers.point(coords[multiFeatureIndex], properties, options);
    case 'LineString':
        if (coordIndex < 0) coordIndex = coords.length + coordIndex;
        return helpers.point(coords[coordIndex], properties, options);
    case 'Polygon':
        if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
        if (coordIndex < 0) coordIndex = coords[geometryIndex].length + coordIndex;
        return helpers.point(coords[geometryIndex][coordIndex], properties, options);
    case 'MultiLineString':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (coordIndex < 0) coordIndex = coords[multiFeatureIndex].length + coordIndex;
        return helpers.point(coords[multiFeatureIndex][coordIndex], properties, options);
    case 'MultiPolygon':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (geometryIndex < 0) geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
        if (coordIndex < 0) coordIndex = coords[multiFeatureIndex][geometryIndex].length - coordIndex;
        return helpers.point(coords[multiFeatureIndex][geometryIndex][coordIndex], properties, options);
    }
    throw new Error('geojson is invalid');
}

exports.coordEach = coordEach;
exports.coordReduce = coordReduce;
exports.propEach = propEach;
exports.propReduce = propReduce;
exports.featureEach = featureEach;
exports.featureReduce = featureReduce;
exports.coordAll = coordAll;
exports.geomEach = geomEach;
exports.geomReduce = geomReduce;
exports.flattenEach = flattenEach;
exports.flattenReduce = flattenReduce;
exports.segmentEach = segmentEach;
exports.segmentReduce = segmentReduce;
exports.lineEach = lineEach;
exports.lineReduce = lineReduce;
exports.findSegment = findSegment;
exports.findPoint = findPoint;

},{"@turf/helpers":15}],12:[function(require,module,exports){
'use strict';

/**
 * Returns a cloned copy of the passed GeoJSON Object, including possible 'Foreign Members'.
 * ~3-5x faster than the common JSON.parse + JSON.stringify combo method.
 *
 * @name clone
 * @param {GeoJSON} geojson GeoJSON Object
 * @returns {GeoJSON} cloned GeoJSON Object
 * @example
 * var line = turf.lineString([[-74, 40], [-78, 42], [-82, 35]], {color: 'red'});
 *
 * var lineCloned = turf.clone(line);
 */
function clone(geojson) {
    if (!geojson) throw new Error('geojson is required');

    switch (geojson.type) {
    case 'Feature':
        return cloneFeature(geojson);
    case 'FeatureCollection':
        return cloneFeatureCollection(geojson);
    case 'Point':
    case 'LineString':
    case 'Polygon':
    case 'MultiPoint':
    case 'MultiLineString':
    case 'MultiPolygon':
    case 'GeometryCollection':
        return cloneGeometry(geojson);
    default:
        throw new Error('unknown GeoJSON type');
    }
}

/**
 * Clone Feature
 *
 * @private
 * @param {Feature<any>} geojson GeoJSON Feature
 * @returns {Feature<any>} cloned Feature
 */
function cloneFeature(geojson) {
    var cloned = {type: 'Feature'};
    // Preserve Foreign Members
    Object.keys(geojson).forEach(function (key) {
        switch (key) {
        case 'type':
        case 'properties':
        case 'geometry':
            return;
        default:
            cloned[key] = geojson[key];
        }
    });
    // Add properties & geometry last
    cloned.properties = cloneProperties(geojson.properties);
    cloned.geometry = cloneGeometry(geojson.geometry);
    return cloned;
}

/**
 * Clone Properties
 *
 * @private
 * @param {Object} properties GeoJSON Properties
 * @returns {Object} cloned Properties
 */
function cloneProperties(properties) {
    var cloned = {};
    if (!properties) return cloned;
    Object.keys(properties).forEach(function (key) {
        var value = properties[key];
        if (typeof value === 'object') {
            if (value === null) {
                // handle null
                cloned[key] = null;
            } else if (value.length) {
                // handle Array
                cloned[key] = value.map(function (item) {
                    return item;
                });
            } else {
                // handle generic Object
                cloned[key] = cloneProperties(value);
            }
        } else cloned[key] = value;
    });
    return cloned;
}

/**
 * Clone Feature Collection
 *
 * @private
 * @param {FeatureCollection<any>} geojson GeoJSON Feature Collection
 * @returns {FeatureCollection<any>} cloned Feature Collection
 */
function cloneFeatureCollection(geojson) {
    var cloned = {type: 'FeatureCollection'};

    // Preserve Foreign Members
    Object.keys(geojson).forEach(function (key) {
        switch (key) {
        case 'type':
        case 'features':
            return;
        default:
            cloned[key] = geojson[key];
        }
    });
    // Add features
    cloned.features = geojson.features.map(function (feature) {
        return cloneFeature(feature);
    });
    return cloned;
}

/**
 * Clone Geometry
 *
 * @private
 * @param {Geometry<any>} geometry GeoJSON Geometry
 * @returns {Geometry<any>} cloned Geometry
 */
function cloneGeometry(geometry) {
    var geom = {type: geometry.type};
    if (geometry.bbox) geom.bbox = geometry.bbox;

    if (geometry.type === 'GeometryCollection') {
        geom.geometries = geometry.geometries.map(function (geom) {
            return cloneGeometry(geom);
        });
        return geom;
    }
    geom.coordinates = deepSlice(geometry.coordinates);
    return geom;
}

/**
 * Deep Slice coordinates
 *
 * @private
 * @param {Coordinates} coords Coordinates
 * @returns {Coordinates} all coordinates sliced
 */
function deepSlice(coords) {
    if (typeof coords[0] !== 'object') { return coords.slice(); }
    return coords.map(function (coord) {
        return deepSlice(coord);
    });
}

module.exports = clone;
module.exports.default = clone;

},{}],13:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
// http://en.wikipedia.org/wiki/Haversine_formula
// http://www.movable-type.co.uk/scripts/latlong.html
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
/**
 * Takes a {@link Point} and calculates the location of a destination point given a distance in
 * degrees, radians, miles, or kilometers; and bearing in degrees.
 * This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
 *
 * @name destination
 * @param {Coord} origin starting point
 * @param {number} distance distance from the origin point
 * @param {number} bearing ranging from -180 to 180
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] miles, kilometers, degrees, or radians
 * @param {Object} [options.properties={}] Translate properties to Point
 * @returns {Feature<Point>} destination point
 * @example
 * var point = turf.point([-75.343, 39.984]);
 * var distance = 50;
 * var bearing = 90;
 * var options = {units: 'miles'};
 *
 * var destination = turf.destination(point, distance, bearing, options);
 *
 * //addToMap
 * var addToMap = [point, destination]
 * destination.properties['marker-color'] = '#f00';
 * point.properties['marker-color'] = '#0f0';
 */
function destination(origin, distance, bearing, options) {
    if (options === void 0) { options = {}; }
    // Handle input
    var coordinates1 = invariant_1.getCoord(origin);
    var longitude1 = helpers_1.degreesToRadians(coordinates1[0]);
    var latitude1 = helpers_1.degreesToRadians(coordinates1[1]);
    var bearingRad = helpers_1.degreesToRadians(bearing);
    var radians = helpers_1.lengthToRadians(distance, options.units);
    // Main
    var latitude2 = Math.asin(Math.sin(latitude1) * Math.cos(radians) +
        Math.cos(latitude1) * Math.sin(radians) * Math.cos(bearingRad));
    var longitude2 = longitude1 + Math.atan2(Math.sin(bearingRad) * Math.sin(radians) * Math.cos(latitude1), Math.cos(radians) - Math.sin(latitude1) * Math.sin(latitude2));
    var lng = helpers_1.radiansToDegrees(longitude2);
    var lat = helpers_1.radiansToDegrees(latitude2);
    return helpers_1.point([lng, lat], options.properties);
}
exports.default = destination;

},{"@turf/helpers":15,"@turf/invariant":17}],14:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var invariant_1 = require("@turf/invariant");
var helpers_1 = require("@turf/helpers");
//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html
/**
 * Calculates the distance between two {@link Point|points} in degrees, radians, miles, or kilometers.
 * This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
 *
 * @name distance
 * @param {Coord} from origin point
 * @param {Coord} to destination point
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
 * @returns {number} distance between the two points
 * @example
 * var from = turf.point([-75.343, 39.984]);
 * var to = turf.point([-75.534, 39.123]);
 * var options = {units: 'miles'};
 *
 * var distance = turf.distance(from, to, options);
 *
 * //addToMap
 * var addToMap = [from, to];
 * from.properties.distance = distance;
 * to.properties.distance = distance;
 */
function distance(from, to, options) {
    if (options === void 0) { options = {}; }
    var coordinates1 = invariant_1.getCoord(from);
    var coordinates2 = invariant_1.getCoord(to);
    var dLat = helpers_1.degreesToRadians((coordinates2[1] - coordinates1[1]));
    var dLon = helpers_1.degreesToRadians((coordinates2[0] - coordinates1[0]));
    var lat1 = helpers_1.degreesToRadians(coordinates1[1]);
    var lat2 = helpers_1.degreesToRadians(coordinates2[1]);
    var a = Math.pow(Math.sin(dLat / 2), 2) +
        Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
    return helpers_1.radiansToLength(2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)), options.units);
}
exports.default = distance;

},{"@turf/helpers":15,"@turf/invariant":17}],15:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * @module helpers
 */
/**
 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
 *
 * @memberof helpers
 * @type {number}
 */
exports.earthRadius = 6371008.8;
/**
 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.factors = {
    centimeters: exports.earthRadius * 100,
    centimetres: exports.earthRadius * 100,
    degrees: exports.earthRadius / 111325,
    feet: exports.earthRadius * 3.28084,
    inches: exports.earthRadius * 39.370,
    kilometers: exports.earthRadius / 1000,
    kilometres: exports.earthRadius / 1000,
    meters: exports.earthRadius,
    metres: exports.earthRadius,
    miles: exports.earthRadius / 1609.344,
    millimeters: exports.earthRadius * 1000,
    millimetres: exports.earthRadius * 1000,
    nauticalmiles: exports.earthRadius / 1852,
    radians: 1,
    yards: exports.earthRadius / 1.0936,
};
/**
 * Units of measurement factors based on 1 meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.unitsFactors = {
    centimeters: 100,
    centimetres: 100,
    degrees: 1 / 111325,
    feet: 3.28084,
    inches: 39.370,
    kilometers: 1 / 1000,
    kilometres: 1 / 1000,
    meters: 1,
    metres: 1,
    miles: 1 / 1609.344,
    millimeters: 1000,
    millimetres: 1000,
    nauticalmiles: 1 / 1852,
    radians: 1 / exports.earthRadius,
    yards: 1 / 1.0936,
};
/**
 * Area of measurement factors based on 1 square meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.areaFactors = {
    acres: 0.000247105,
    centimeters: 10000,
    centimetres: 10000,
    feet: 10.763910417,
    inches: 1550.003100006,
    kilometers: 0.000001,
    kilometres: 0.000001,
    meters: 1,
    metres: 1,
    miles: 3.86e-7,
    millimeters: 1000000,
    millimetres: 1000000,
    yards: 1.195990046,
};
/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature} a GeoJSON Feature
 * @example
 * var geometry = {
 *   "type": "Point",
 *   "coordinates": [110, 50]
 * };
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geom, properties, options) {
    if (options === void 0) { options = {}; }
    var feat = { type: "Feature" };
    if (options.id === 0 || options.id) {
        feat.id = options.id;
    }
    if (options.bbox) {
        feat.bbox = options.bbox;
    }
    feat.properties = properties || {};
    feat.geometry = geom;
    return feat;
}
exports.feature = feature;
/**
 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
 * For GeometryCollection type use `helpers.geometryCollection`
 *
 * @name geometry
 * @param {string} type Geometry Type
 * @param {Array<any>} coordinates Coordinates
 * @param {Object} [options={}] Optional Parameters
 * @returns {Geometry} a GeoJSON Geometry
 * @example
 * var type = "Point";
 * var coordinates = [110, 50];
 * var geometry = turf.geometry(type, coordinates);
 * // => geometry
 */
function geometry(type, coordinates, options) {
    if (options === void 0) { options = {}; }
    switch (type) {
        case "Point": return point(coordinates).geometry;
        case "LineString": return lineString(coordinates).geometry;
        case "Polygon": return polygon(coordinates).geometry;
        case "MultiPoint": return multiPoint(coordinates).geometry;
        case "MultiLineString": return multiLineString(coordinates).geometry;
        case "MultiPolygon": return multiPolygon(coordinates).geometry;
        default: throw new Error(type + " is invalid");
    }
}
exports.geometry = geometry;
/**
 * Creates a {@link Point} {@link Feature} from a Position.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Point>} a Point feature
 * @example
 * var point = turf.point([-75.343, 39.984]);
 *
 * //=point
 */
function point(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "Point",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.point = point;
/**
 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
 *
 * @name points
 * @param {Array<Array<number>>} coordinates an array of Points
 * @param {Object} [properties={}] Translate these properties to each Feature
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Point>} Point Feature
 * @example
 * var points = turf.points([
 *   [-75, 39],
 *   [-80, 45],
 *   [-78, 50]
 * ]);
 *
 * //=points
 */
function points(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return point(coords, properties);
    }), options);
}
exports.points = points;
/**
 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Polygon>} Polygon Feature
 * @example
 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
 *
 * //=polygon
 */
function polygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    for (var _i = 0, coordinates_1 = coordinates; _i < coordinates_1.length; _i++) {
        var ring = coordinates_1[_i];
        if (ring.length < 4) {
            throw new Error("Each LinearRing of a Polygon must have 4 or more Positions.");
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            // Check if first point of Polygon contains two numbers
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error("First and last Position are not equivalent.");
            }
        }
    }
    var geom = {
        type: "Polygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.polygon = polygon;
/**
 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
 *
 * @name polygons
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
 * @example
 * var polygons = turf.polygons([
 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
 * ]);
 *
 * //=polygons
 */
function polygons(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return polygon(coords, properties);
    }), options);
}
exports.polygons = polygons;
/**
 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<LineString>} LineString Feature
 * @example
 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
 *
 * //=linestring1
 * //=linestring2
 */
function lineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    if (coordinates.length < 2) {
        throw new Error("coordinates must be an array of two or more positions");
    }
    var geom = {
        type: "LineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.lineString = lineString;
/**
 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
 *
 * @name lineStrings
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
 * @example
 * var linestrings = turf.lineStrings([
 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
 * ]);
 *
 * //=linestrings
 */
function lineStrings(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return lineString(coords, properties);
    }), options);
}
exports.lineStrings = lineStrings;
/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {FeatureCollection} FeatureCollection of Features
 * @example
 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
 *
 * var collection = turf.featureCollection([
 *   locationA,
 *   locationB,
 *   locationC
 * ]);
 *
 * //=collection
 */
function featureCollection(features, options) {
    if (options === void 0) { options = {}; }
    var fc = { type: "FeatureCollection" };
    if (options.id) {
        fc.id = options.id;
    }
    if (options.bbox) {
        fc.bbox = options.bbox;
    }
    fc.features = features;
    return fc;
}
exports.featureCollection = featureCollection;
/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 */
function multiLineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiLineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiLineString = multiLineString;
/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 */
function multiPoint(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPoint",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPoint = multiPoint;
/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
function multiPolygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPolygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPolygon = multiPolygon;
/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = turf.geometry("Point", [100, 0]);
 * var line = turf.geometry("LineString", [[101, 0], [102, 1]]);
 * var collection = turf.geometryCollection([pt, line]);
 *
 * // => collection
 */
function geometryCollection(geometries, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "GeometryCollection",
        geometries: geometries,
    };
    return feature(geom, properties, options);
}
exports.geometryCollection = geometryCollection;
/**
 * Round number to precision
 *
 * @param {number} num Number
 * @param {number} [precision=0] Precision
 * @returns {number} rounded number
 * @example
 * turf.round(120.4321)
 * //=120
 *
 * turf.round(120.4321, 2)
 * //=120.43
 */
function round(num, precision) {
    if (precision === void 0) { precision = 0; }
    if (precision && !(precision >= 0)) {
        throw new Error("precision must be a positive number");
    }
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(num * multiplier) / multiplier;
}
exports.round = round;
/**
 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name radiansToLength
 * @param {number} radians in radians across the sphere
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} distance
 */
function radiansToLength(radians, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return radians * factor;
}
exports.radiansToLength = radiansToLength;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name lengthToRadians
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} radians
 */
function lengthToRadians(distance, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return distance / factor;
}
exports.lengthToRadians = lengthToRadians;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
 *
 * @name lengthToDegrees
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} degrees
 */
function lengthToDegrees(distance, units) {
    return radiansToDegrees(lengthToRadians(distance, units));
}
exports.lengthToDegrees = lengthToDegrees;
/**
 * Converts any bearing angle from the north line direction (positive clockwise)
 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
 *
 * @name bearingToAzimuth
 * @param {number} bearing angle, between -180 and +180 degrees
 * @returns {number} angle between 0 and 360 degrees
 */
function bearingToAzimuth(bearing) {
    var angle = bearing % 360;
    if (angle < 0) {
        angle += 360;
    }
    return angle;
}
exports.bearingToAzimuth = bearingToAzimuth;
/**
 * Converts an angle in radians to degrees
 *
 * @name radiansToDegrees
 * @param {number} radians angle in radians
 * @returns {number} degrees between 0 and 360 degrees
 */
function radiansToDegrees(radians) {
    var degrees = radians % (2 * Math.PI);
    return degrees * 180 / Math.PI;
}
exports.radiansToDegrees = radiansToDegrees;
/**
 * Converts an angle in degrees to radians
 *
 * @name degreesToRadians
 * @param {number} degrees angle between 0 and 360 degrees
 * @returns {number} angle in radians
 */
function degreesToRadians(degrees) {
    var radians = degrees % 360;
    return radians * Math.PI / 180;
}
exports.degreesToRadians = degreesToRadians;
/**
 * Converts a length to the requested unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @param {number} length to be converted
 * @param {Units} [originalUnit="kilometers"] of the length
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted length
 */
function convertLength(length, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "kilometers"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(length >= 0)) {
        throw new Error("length must be a positive number");
    }
    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
}
exports.convertLength = convertLength;
/**
 * Converts a area to the requested unit.
 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches
 * @param {number} area to be converted
 * @param {Units} [originalUnit="meters"] of the distance
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted distance
 */
function convertArea(area, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "meters"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(area >= 0)) {
        throw new Error("area must be a positive number");
    }
    var startFactor = exports.areaFactors[originalUnit];
    if (!startFactor) {
        throw new Error("invalid original units");
    }
    var finalFactor = exports.areaFactors[finalUnit];
    if (!finalFactor) {
        throw new Error("invalid final units");
    }
    return (area / startFactor) * finalFactor;
}
exports.convertArea = convertArea;
/**
 * isNumber
 *
 * @param {*} num Number to validate
 * @returns {boolean} true/false
 * @example
 * turf.isNumber(123)
 * //=true
 * turf.isNumber('foo')
 * //=false
 */
function isNumber(num) {
    return !isNaN(num) && num !== null && !Array.isArray(num) && !/^\s*$/.test(num);
}
exports.isNumber = isNumber;
/**
 * isObject
 *
 * @param {*} input variable to validate
 * @returns {boolean} true/false
 * @example
 * turf.isObject({elevation: 10})
 * //=true
 * turf.isObject('foo')
 * //=false
 */
function isObject(input) {
    return (!!input) && (input.constructor === Object);
}
exports.isObject = isObject;
/**
 * Validate BBox
 *
 * @private
 * @param {Array<number>} bbox BBox to validate
 * @returns {void}
 * @throws Error if BBox is not valid
 * @example
 * validateBBox([-180, -40, 110, 50])
 * //=OK
 * validateBBox([-180, -40])
 * //=Error
 * validateBBox('Foo')
 * //=Error
 * validateBBox(5)
 * //=Error
 * validateBBox(null)
 * //=Error
 * validateBBox(undefined)
 * //=Error
 */
function validateBBox(bbox) {
    if (!bbox) {
        throw new Error("bbox is required");
    }
    if (!Array.isArray(bbox)) {
        throw new Error("bbox must be an Array");
    }
    if (bbox.length !== 4 && bbox.length !== 6) {
        throw new Error("bbox must be an Array of 4 or 6 numbers");
    }
    bbox.forEach(function (num) {
        if (!isNumber(num)) {
            throw new Error("bbox must only contain numbers");
        }
    });
}
exports.validateBBox = validateBBox;
/**
 * Validate Id
 *
 * @private
 * @param {string|number} id Id to validate
 * @returns {void}
 * @throws Error if Id is not valid
 * @example
 * validateId([-180, -40, 110, 50])
 * //=Error
 * validateId([-180, -40])
 * //=Error
 * validateId('Foo')
 * //=OK
 * validateId(5)
 * //=OK
 * validateId(null)
 * //=Error
 * validateId(undefined)
 * //=Error
 */
function validateId(id) {
    if (!id) {
        throw new Error("id is required");
    }
    if (["string", "number"].indexOf(typeof id) === -1) {
        throw new Error("id must be a number or a string");
    }
}
exports.validateId = validateId;
// Deprecated methods
function radians2degrees() {
    throw new Error("method has been renamed to `radiansToDegrees`");
}
exports.radians2degrees = radians2degrees;
function degrees2radians() {
    throw new Error("method has been renamed to `degreesToRadians`");
}
exports.degrees2radians = degrees2radians;
function distanceToDegrees() {
    throw new Error("method has been renamed to `lengthToDegrees`");
}
exports.distanceToDegrees = distanceToDegrees;
function distanceToRadians() {
    throw new Error("method has been renamed to `lengthToRadians`");
}
exports.distanceToRadians = distanceToRadians;
function radiansToDistance() {
    throw new Error("method has been renamed to `radiansToLength`");
}
exports.radiansToDistance = radiansToDistance;
function bearingToAngle() {
    throw new Error("method has been renamed to `bearingToAzimuth`");
}
exports.bearingToAngle = bearingToAngle;
function convertDistance() {
    throw new Error("method has been renamed to `convertLength`");
}
exports.convertDistance = convertDistance;

},{}],16:[function(require,module,exports){
"use strict";
var __importStar = (this && this.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (Object.hasOwnProperty.call(mod, k)) result[k] = mod[k];
    result["default"] = mod;
    return result;
};
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
var martinez = __importStar(require("martinez-polygon-clipping"));
/**
 * Takes two {@link Polygon|polygon} or {@link MultiPolygon|multi-polygon} geometries and
 * finds their polygonal intersection. If they don't intersect, returns null.
 *
 * @name intersect
 * @param {Feature<Polygon | MultiPolygon>} poly1 the first polygon or multipolygon
 * @param {Feature<Polygon | MultiPolygon>} poly2 the second polygon or multipolygon
 * @param {Object} [options={}] Optional Parameters
 * @param {Object} [options.properties={}] Translate GeoJSON Properties to Feature
 * @returns {Feature|null} returns a feature representing the area they share (either a {@link Polygon} or
 * {@link MultiPolygon}). If they do not share any area, returns `null`.
 * @example
 * var poly1 = turf.polygon([[
 *   [-122.801742, 45.48565],
 *   [-122.801742, 45.60491],
 *   [-122.584762, 45.60491],
 *   [-122.584762, 45.48565],
 *   [-122.801742, 45.48565]
 * ]]);
 *
 * var poly2 = turf.polygon([[
 *   [-122.520217, 45.535693],
 *   [-122.64038, 45.553967],
 *   [-122.720031, 45.526554],
 *   [-122.669906, 45.507309],
 *   [-122.723464, 45.446643],
 *   [-122.532577, 45.408574],
 *   [-122.487258, 45.477466],
 *   [-122.520217, 45.535693]
 * ]]);
 *
 * var intersection = turf.intersect(poly1, poly2);
 *
 * //addToMap
 * var addToMap = [poly1, poly2, intersection];
 */
function intersect(poly1, poly2, options) {
    if (options === void 0) { options = {}; }
    var geom1 = invariant_1.getGeom(poly1);
    var geom2 = invariant_1.getGeom(poly2);
    if (geom1.type === "Polygon" && geom2.type === "Polygon") {
        var intersection = martinez.intersection(geom1.coordinates, geom2.coordinates);
        if (intersection === null || intersection.length === 0) {
            return null;
        }
        if (intersection.length === 1) {
            var start = intersection[0][0][0];
            var end = intersection[0][0][intersection[0][0].length - 1];
            if (start[0] === end[0] && start[1] === end[1]) {
                return helpers_1.polygon(intersection[0], options.properties);
            }
            return null;
        }
        return helpers_1.multiPolygon(intersection, options.properties);
    }
    else if (geom1.type === "MultiPolygon") {
        var resultCoords = [];
        // iterate through the polygon and run intersect with each part, adding to the resultCoords.
        for (var _i = 0, _a = geom1.coordinates; _i < _a.length; _i++) {
            var coords = _a[_i];
            var subGeom = invariant_1.getGeom(helpers_1.polygon(coords));
            var subIntersection = intersect(subGeom, geom2);
            if (subIntersection) {
                var subIntGeom = invariant_1.getGeom(subIntersection);
                if (subIntGeom.type === "Polygon") {
                    resultCoords.push(subIntGeom.coordinates);
                }
                else if (subIntGeom.type === "MultiPolygon") {
                    resultCoords = resultCoords.concat(subIntGeom.coordinates);
                }
                else {
                    throw new Error("intersection is invalid");
                }
            }
        }
        // Make a polygon with the result
        if (resultCoords.length === 0) {
            return null;
        }
        if (resultCoords.length === 1) {
            return helpers_1.polygon(resultCoords[0], options.properties);
        }
        else {
            return helpers_1.multiPolygon(resultCoords, options.properties);
        }
    }
    else if (geom2.type === "MultiPolygon") {
        // geom1 is a polygon and geom2 a multiPolygon,
        // put the multiPolygon first and fallback to the previous case.
        return intersect(geom2, geom1);
    }
    else {
        // handle invalid geometry types
        throw new Error("poly1 and poly2 must be either polygons or multiPolygons");
    }
}
exports.default = intersect;

},{"@turf/helpers":15,"@turf/invariant":17,"martinez-polygon-clipping":53}],17:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("@turf/helpers");
/**
 * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
 *
 * @name getCoord
 * @param {Array<number>|Geometry<Point>|Feature<Point>} coord GeoJSON Point or an Array of numbers
 * @returns {Array<number>} coordinates
 * @example
 * var pt = turf.point([10, 10]);
 *
 * var coord = turf.getCoord(pt);
 * //= [10, 10]
 */
function getCoord(coord) {
    if (!coord) {
        throw new Error("coord is required");
    }
    if (!Array.isArray(coord)) {
        if (coord.type === "Feature" && coord.geometry !== null && coord.geometry.type === "Point") {
            return coord.geometry.coordinates;
        }
        if (coord.type === "Point") {
            return coord.coordinates;
        }
    }
    if (Array.isArray(coord) && coord.length >= 2 && !Array.isArray(coord[0]) && !Array.isArray(coord[1])) {
        return coord;
    }
    throw new Error("coord must be GeoJSON Point or an Array of numbers");
}
exports.getCoord = getCoord;
/**
 * Unwrap coordinates from a Feature, Geometry Object or an Array
 *
 * @name getCoords
 * @param {Array<any>|Geometry|Feature} coords Feature, Geometry Object or an Array
 * @returns {Array<any>} coordinates
 * @example
 * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
 *
 * var coords = turf.getCoords(poly);
 * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
 */
function getCoords(coords) {
    if (Array.isArray(coords)) {
        return coords;
    }
    // Feature
    if (coords.type === "Feature") {
        if (coords.geometry !== null) {
            return coords.geometry.coordinates;
        }
    }
    else {
        // Geometry
        if (coords.coordinates) {
            return coords.coordinates;
        }
    }
    throw new Error("coords must be GeoJSON Feature, Geometry Object or an Array");
}
exports.getCoords = getCoords;
/**
 * Checks if coordinates contains a number
 *
 * @name containsNumber
 * @param {Array<any>} coordinates GeoJSON Coordinates
 * @returns {boolean} true if Array contains a number
 */
function containsNumber(coordinates) {
    if (coordinates.length > 1 && helpers_1.isNumber(coordinates[0]) && helpers_1.isNumber(coordinates[1])) {
        return true;
    }
    if (Array.isArray(coordinates[0]) && coordinates[0].length) {
        return containsNumber(coordinates[0]);
    }
    throw new Error("coordinates must only contain numbers");
}
exports.containsNumber = containsNumber;
/**
 * Enforce expectations about types of GeoJSON objects for Turf.
 *
 * @name geojsonType
 * @param {GeoJSON} value any GeoJSON object
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function geojsonType(value, type, name) {
    if (!type || !name) {
        throw new Error("type and name required");
    }
    if (!value || value.type !== type) {
        throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + value.type);
    }
}
exports.geojsonType = geojsonType;
/**
 * Enforce expectations about types of {@link Feature} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name featureOf
 * @param {Feature} feature a feature with an expected geometry type
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} error if value is not the expected type.
 */
function featureOf(feature, type, name) {
    if (!feature) {
        throw new Error("No feature passed");
    }
    if (!name) {
        throw new Error(".featureOf() requires a name");
    }
    if (!feature || feature.type !== "Feature" || !feature.geometry) {
        throw new Error("Invalid input to " + name + ", Feature with geometry required");
    }
    if (!feature.geometry || feature.geometry.type !== type) {
        throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + feature.geometry.type);
    }
}
exports.featureOf = featureOf;
/**
 * Enforce expectations about types of {@link FeatureCollection} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name collectionOf
 * @param {FeatureCollection} featureCollection a FeatureCollection for which features will be judged
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function collectionOf(featureCollection, type, name) {
    if (!featureCollection) {
        throw new Error("No featureCollection passed");
    }
    if (!name) {
        throw new Error(".collectionOf() requires a name");
    }
    if (!featureCollection || featureCollection.type !== "FeatureCollection") {
        throw new Error("Invalid input to " + name + ", FeatureCollection required");
    }
    for (var _i = 0, _a = featureCollection.features; _i < _a.length; _i++) {
        var feature = _a[_i];
        if (!feature || feature.type !== "Feature" || !feature.geometry) {
            throw new Error("Invalid input to " + name + ", Feature with geometry required");
        }
        if (!feature.geometry || feature.geometry.type !== type) {
            throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + feature.geometry.type);
        }
    }
}
exports.collectionOf = collectionOf;
/**
 * Get Geometry from Feature or Geometry Object
 *
 * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
 * @returns {Geometry|null} GeoJSON Geometry Object
 * @throws {Error} if geojson is not a Feature or Geometry Object
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getGeom(point)
 * //={"type": "Point", "coordinates": [110, 40]}
 */
function getGeom(geojson) {
    if (geojson.type === "Feature") {
        return geojson.geometry;
    }
    return geojson;
}
exports.getGeom = getGeom;
/**
 * Get GeoJSON object's type, Geometry type is prioritize.
 *
 * @param {GeoJSON} geojson GeoJSON object
 * @param {string} [name="geojson"] name of the variable to display in error message
 * @returns {string} GeoJSON type
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getType(point)
 * //="Point"
 */
function getType(geojson, name) {
    if (geojson.type === "FeatureCollection") {
        return "FeatureCollection";
    }
    if (geojson.type === "GeometryCollection") {
        return "GeometryCollection";
    }
    if (geojson.type === "Feature" && geojson.geometry !== null) {
        return geojson.geometry.type;
    }
    return geojson.type;
}
exports.getType = getType;

},{"@turf/helpers":15}],18:[function(require,module,exports){
"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
}
Object.defineProperty(exports, "__esModule", { value: true });
var distance_1 = __importDefault(require("@turf/distance"));
var meta_1 = require("@turf/meta");
/**
 * Takes a {@link GeoJSON} and measures its length in the specified units, {@link (Multi)Point}'s distance are ignored.
 *
 * @name length
 * @param {Feature<LineString|MultiLineString>} geojson GeoJSON to measure
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units=kilometers] can be degrees, radians, miles, or kilometers
 * @returns {number} length of GeoJSON
 * @example
 * var line = turf.lineString([[115, -32], [131, -22], [143, -25], [150, -34]]);
 * var length = turf.length(line, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [line];
 * line.properties.distance = length;
 */
function length(geojson, options) {
    if (options === void 0) { options = {}; }
    // Calculate distance from 2-vertex line segments
    return meta_1.segmentReduce(geojson, function (previousValue, segment) {
        var coords = segment.geometry.coordinates;
        return previousValue + distance_1.default(coords[0], coords[1], options);
    }, 0);
}
exports.default = length;

},{"@turf/distance":14,"@turf/meta":19}],19:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],20:[function(require,module,exports){
"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
}
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
var line_segment_1 = __importDefault(require("@turf/line-segment"));
var meta_1 = require("@turf/meta");
var geojson_rbush_1 = __importDefault(require("geojson-rbush"));
/**
 * Takes any LineString or Polygon GeoJSON and returns the intersecting point(s).
 *
 * @name lineIntersect
 * @param {GeoJSON} line1 any LineString or Polygon
 * @param {GeoJSON} line2 any LineString or Polygon
 * @returns {FeatureCollection<Point>} point(s) that intersect both
 * @example
 * var line1 = turf.lineString([[126, -11], [129, -21]]);
 * var line2 = turf.lineString([[123, -18], [131, -14]]);
 * var intersects = turf.lineIntersect(line1, line2);
 *
 * //addToMap
 * var addToMap = [line1, line2, intersects]
 */
function lineIntersect(line1, line2) {
    var unique = {};
    var results = [];
    // First, normalize geometries to features
    // Then, handle simple 2-vertex segments
    if (line1.type === "LineString") {
        line1 = helpers_1.feature(line1);
    }
    if (line2.type === "LineString") {
        line2 = helpers_1.feature(line2);
    }
    if (line1.type === "Feature" &&
        line2.type === "Feature" &&
        line1.geometry !== null &&
        line2.geometry !== null &&
        line1.geometry.type === "LineString" &&
        line2.geometry.type === "LineString" &&
        line1.geometry.coordinates.length === 2 &&
        line2.geometry.coordinates.length === 2) {
        var intersect = intersects(line1, line2);
        if (intersect) {
            results.push(intersect);
        }
        return helpers_1.featureCollection(results);
    }
    // Handles complex GeoJSON Geometries
    var tree = geojson_rbush_1.default();
    tree.load(line_segment_1.default(line2));
    meta_1.featureEach(line_segment_1.default(line1), function (segment) {
        meta_1.featureEach(tree.search(segment), function (match) {
            var intersect = intersects(segment, match);
            if (intersect) {
                // prevent duplicate points https://github.com/Turfjs/turf/issues/688
                var key = invariant_1.getCoords(intersect).join(",");
                if (!unique[key]) {
                    unique[key] = true;
                    results.push(intersect);
                }
            }
        });
    });
    return helpers_1.featureCollection(results);
}
/**
 * Find a point that intersects LineStrings with two coordinates each
 *
 * @private
 * @param {Feature<LineString>} line1 GeoJSON LineString (Must only contain 2 coordinates)
 * @param {Feature<LineString>} line2 GeoJSON LineString (Must only contain 2 coordinates)
 * @returns {Feature<Point>} intersecting GeoJSON Point
 */
function intersects(line1, line2) {
    var coords1 = invariant_1.getCoords(line1);
    var coords2 = invariant_1.getCoords(line2);
    if (coords1.length !== 2) {
        throw new Error("<intersects> line1 must only contain 2 coordinates");
    }
    if (coords2.length !== 2) {
        throw new Error("<intersects> line2 must only contain 2 coordinates");
    }
    var x1 = coords1[0][0];
    var y1 = coords1[0][1];
    var x2 = coords1[1][0];
    var y2 = coords1[1][1];
    var x3 = coords2[0][0];
    var y3 = coords2[0][1];
    var x4 = coords2[1][0];
    var y4 = coords2[1][1];
    var denom = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
    var numeA = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
    var numeB = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));
    if (denom === 0) {
        if (numeA === 0 && numeB === 0) {
            return null;
        }
        return null;
    }
    var uA = numeA / denom;
    var uB = numeB / denom;
    if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
        var x = x1 + (uA * (x2 - x1));
        var y = y1 + (uA * (y2 - y1));
        return helpers_1.point([x, y]);
    }
    return null;
}
exports.default = lineIntersect;

},{"@turf/helpers":15,"@turf/invariant":17,"@turf/line-segment":22,"@turf/meta":21,"geojson-rbush":51}],21:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],22:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
var meta_1 = require("@turf/meta");
/**
 * Creates a {@link FeatureCollection} of 2-vertex {@link LineString} segments from a
 * {@link LineString|(Multi)LineString} or {@link Polygon|(Multi)Polygon}.
 *
 * @name lineSegment
 * @param {GeoJSON} geojson GeoJSON Polygon or LineString
 * @returns {FeatureCollection<LineString>} 2-vertex line segments
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 * var segments = turf.lineSegment(polygon);
 *
 * //addToMap
 * var addToMap = [polygon, segments]
 */
function lineSegment(geojson) {
    if (!geojson) {
        throw new Error("geojson is required");
    }
    var results = [];
    meta_1.flattenEach(geojson, function (feature) {
        lineSegmentFeature(feature, results);
    });
    return helpers_1.featureCollection(results);
}
/**
 * Line Segment
 *
 * @private
 * @param {Feature<LineString|Polygon>} geojson Line or polygon feature
 * @param {Array} results push to results
 * @returns {void}
 */
function lineSegmentFeature(geojson, results) {
    var coords = [];
    var geometry = geojson.geometry;
    if (geometry !== null) {
        switch (geometry.type) {
            case "Polygon":
                coords = invariant_1.getCoords(geometry);
                break;
            case "LineString":
                coords = [invariant_1.getCoords(geometry)];
        }
        coords.forEach(function (coord) {
            var segments = createSegments(coord, geojson.properties);
            segments.forEach(function (segment) {
                segment.id = results.length;
                results.push(segment);
            });
        });
    }
}
/**
 * Create Segments from LineString coordinates
 *
 * @private
 * @param {Array<Array<number>>} coords LineString coordinates
 * @param {*} properties GeoJSON properties
 * @returns {Array<Feature<LineString>>} line segments
 */
function createSegments(coords, properties) {
    var segments = [];
    coords.reduce(function (previousCoords, currentCoords) {
        var segment = helpers_1.lineString([previousCoords, currentCoords], properties);
        segment.bbox = bbox(previousCoords, currentCoords);
        segments.push(segment);
        return currentCoords;
    });
    return segments;
}
/**
 * Create BBox between two coordinates (faster than @turf/bbox)
 *
 * @private
 * @param {Array<number>} coords1 Point coordinate
 * @param {Array<number>} coords2 Point coordinate
 * @returns {BBox} [west, south, east, north]
 */
function bbox(coords1, coords2) {
    var x1 = coords1[0];
    var y1 = coords1[1];
    var x2 = coords2[0];
    var y2 = coords2[1];
    var west = (x1 < x2) ? x1 : x2;
    var south = (y1 < y2) ? y1 : y2;
    var east = (x1 > x2) ? x1 : x2;
    var north = (y1 > y2) ? y1 : y2;
    return [west, south, east, north];
}
exports.default = lineSegment;

},{"@turf/helpers":15,"@turf/invariant":17,"@turf/meta":23}],23:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],24:[function(require,module,exports){
'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var invariant = require('@turf/invariant');
var helpers = require('@turf/helpers');
var nearestPointOnLine = _interopDefault(require('@turf/nearest-point-on-line'));

/**
 * Takes a {@link LineString|line}, a start {@link Point}, and a stop point
 * and returns a subsection of the line in-between those points.
 * The start & stop points don't need to fall exactly on the line.
 *
 * This can be useful for extracting only the part of a route between waypoints.
 *
 * @name lineSlice
 * @param {Coord} startPt starting point
 * @param {Coord} stopPt stopping point
 * @param {Feature<LineString>|LineString} line line to slice
 * @returns {Feature<LineString>} sliced line
 * @example
 * var line = turf.lineString([
 *     [-77.031669, 38.878605],
 *     [-77.029609, 38.881946],
 *     [-77.020339, 38.884084],
 *     [-77.025661, 38.885821],
 *     [-77.021884, 38.889563],
 *     [-77.019824, 38.892368]
 * ]);
 * var start = turf.point([-77.029609, 38.881946]);
 * var stop = turf.point([-77.021884, 38.889563]);
 *
 * var sliced = turf.lineSlice(start, stop, line);
 *
 * //addToMap
 * var addToMap = [start, stop, line]
 */
function lineSlice(startPt, stopPt, line) {
    // Validation
    var coords = invariant.getCoords(line);
    if (invariant.getType(line) !== 'LineString') throw new Error('line must be a LineString');

    var startVertex = nearestPointOnLine(line, startPt);
    var stopVertex = nearestPointOnLine(line, stopPt);
    var ends;
    if (startVertex.properties.index <= stopVertex.properties.index) {
        ends = [startVertex, stopVertex];
    } else {
        ends = [stopVertex, startVertex];
    }
    var clipCoords = [ends[0].geometry.coordinates];
    for (var i = ends[0].properties.index + 1; i < ends[1].properties.index + 1; i++) {
        clipCoords.push(coords[i]);
    }
    clipCoords.push(ends[1].geometry.coordinates);
    return helpers.lineString(clipCoords, line.properties);
}

module.exports = lineSlice;
module.exports.default = lineSlice;

},{"@turf/helpers":28,"@turf/invariant":29,"@turf/nearest-point-on-line":32}],25:[function(require,module,exports){
'use strict';

var invariant = require('@turf/invariant');
var helpers = require('@turf/helpers');

//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html

/**
 * Takes two {@link Point|points} and finds the geographic bearing between them,
 * i.e. the angle measured in degrees from the north line (0 degrees)
 *
 * @name bearing
 * @param {Coord} start starting Point
 * @param {Coord} end ending Point
 * @param {Object} [options={}] Optional parameters
 * @param {boolean} [options.final=false] calculates the final bearing if true
 * @returns {number} bearing in decimal degrees, between -180 and 180 degrees (positive clockwise)
 * @example
 * var point1 = turf.point([-75.343, 39.984]);
 * var point2 = turf.point([-75.534, 39.123]);
 *
 * var bearing = turf.bearing(point1, point2);
 *
 * //addToMap
 * var addToMap = [point1, point2]
 * point1.properties['marker-color'] = '#f00'
 * point2.properties['marker-color'] = '#0f0'
 * point1.properties.bearing = bearing
 */
function bearing(start, end, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var final = options.final;

    // Reverse calculation
    if (final === true) return calculateFinalBearing(start, end);

    var coordinates1 = invariant.getCoord(start);
    var coordinates2 = invariant.getCoord(end);

    var lon1 = helpers.degreesToRadians(coordinates1[0]);
    var lon2 = helpers.degreesToRadians(coordinates2[0]);
    var lat1 = helpers.degreesToRadians(coordinates1[1]);
    var lat2 = helpers.degreesToRadians(coordinates2[1]);
    var a = Math.sin(lon2 - lon1) * Math.cos(lat2);
    var b = Math.cos(lat1) * Math.sin(lat2) -
        Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);

    return helpers.radiansToDegrees(Math.atan2(a, b));
}

/**
 * Calculates Final Bearing
 *
 * @private
 * @param {Coord} start starting Point
 * @param {Coord} end ending Point
 * @returns {number} bearing
 */
function calculateFinalBearing(start, end) {
    // Swap start & end
    var bear = bearing(end, start);
    bear = (bear + 180) % 360;
    return bear;
}

module.exports = bearing;
module.exports.default = bearing;

},{"@turf/helpers":28,"@turf/invariant":29}],26:[function(require,module,exports){
'use strict';

var invariant = require('@turf/invariant');
var helpers = require('@turf/helpers');

//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html
/**
 * Takes a {@link Point} and calculates the location of a destination point given a distance in degrees, radians, miles, or kilometers; and bearing in degrees. This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
 *
 * @name destination
 * @param {Coord} origin starting point
 * @param {number} distance distance from the origin point
 * @param {number} bearing ranging from -180 to 180
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] miles, kilometers, degrees, or radians
 * @param {Object} [options.properties={}] Translate properties to Point
 * @returns {Feature<Point>} destination point
 * @example
 * var point = turf.point([-75.343, 39.984]);
 * var distance = 50;
 * var bearing = 90;
 * var options = {units: 'miles'};
 *
 * var destination = turf.destination(point, distance, bearing, options);
 *
 * //addToMap
 * var addToMap = [point, destination]
 * destination.properties['marker-color'] = '#f00';
 * point.properties['marker-color'] = '#0f0';
 */
function destination(origin, distance, bearing, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var units = options.units;
    var properties = options.properties;

    // Handle input
    var coordinates1 = invariant.getCoord(origin);
    var longitude1 = helpers.degreesToRadians(coordinates1[0]);
    var latitude1 = helpers.degreesToRadians(coordinates1[1]);
    var bearing_rad = helpers.degreesToRadians(bearing);
    var radians = helpers.lengthToRadians(distance, units);

    // Main
    var latitude2 = Math.asin(Math.sin(latitude1) * Math.cos(radians) +
        Math.cos(latitude1) * Math.sin(radians) * Math.cos(bearing_rad));
    var longitude2 = longitude1 + Math.atan2(Math.sin(bearing_rad) * Math.sin(radians) * Math.cos(latitude1),
        Math.cos(radians) - Math.sin(latitude1) * Math.sin(latitude2));
    var lng = helpers.radiansToDegrees(longitude2);
    var lat = helpers.radiansToDegrees(latitude2);

    return helpers.point([lng, lat], properties);
}

module.exports = destination;
module.exports.default = destination;

},{"@turf/helpers":28,"@turf/invariant":29}],27:[function(require,module,exports){
'use strict';

var invariant = require('@turf/invariant');
var helpers = require('@turf/helpers');

//http://en.wikipedia.org/wiki/Haversine_formula
//http://www.movable-type.co.uk/scripts/latlong.html

/**
 * Calculates the distance between two {@link Point|points} in degrees, radians,
 * miles, or kilometers. This uses the
 * [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula)
 * to account for global curvature.
 *
 * @name distance
 * @param {Coord} from origin point
 * @param {Coord} to destination point
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
 * @returns {number} distance between the two points
 * @example
 * var from = turf.point([-75.343, 39.984]);
 * var to = turf.point([-75.534, 39.123]);
 * var options = {units: 'miles'};
 *
 * var distance = turf.distance(from, to, options);
 *
 * //addToMap
 * var addToMap = [from, to];
 * from.properties.distance = distance;
 * to.properties.distance = distance;
 */
function distance(from, to, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var units = options.units;

    var coordinates1 = invariant.getCoord(from);
    var coordinates2 = invariant.getCoord(to);
    var dLat = helpers.degreesToRadians((coordinates2[1] - coordinates1[1]));
    var dLon = helpers.degreesToRadians((coordinates2[0] - coordinates1[0]));
    var lat1 = helpers.degreesToRadians(coordinates1[1]);
    var lat2 = helpers.degreesToRadians(coordinates2[1]);

    var a = Math.pow(Math.sin(dLat / 2), 2) +
          Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);

    return helpers.radiansToLength(2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)), units);
}

module.exports = distance;
module.exports.default = distance;

},{"@turf/helpers":28,"@turf/invariant":29}],28:[function(require,module,exports){
arguments[4][7][0].apply(exports,arguments)
},{"dup":7}],29:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

var helpers = require('@turf/helpers');

/**
 * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
 *
 * @name getCoord
 * @param {Array<number>|Geometry<Point>|Feature<Point>} coord GeoJSON Point or an Array of numbers
 * @returns {Array<number>} coordinates
 * @example
 * var pt = turf.point([10, 10]);
 *
 * var coord = turf.getCoord(pt);
 * //= [10, 10]
 */
function getCoord(coord) {
    if (!coord) throw new Error('coord is required');
    if (coord.type === 'Feature' && coord.geometry !== null && coord.geometry.type === 'Point') return coord.geometry.coordinates;
    if (coord.type === 'Point') return coord.coordinates;
    if (Array.isArray(coord) && coord.length >= 2 && coord[0].length === undefined && coord[1].length === undefined) return coord;

    throw new Error('coord must be GeoJSON Point or an Array of numbers');
}

/**
 * Unwrap coordinates from a Feature, Geometry Object or an Array
 *
 * @name getCoords
 * @param {Array<any>|Geometry|Feature} coords Feature, Geometry Object or an Array
 * @returns {Array<any>} coordinates
 * @example
 * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
 *
 * var coords = turf.getCoords(poly);
 * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
 */
function getCoords(coords) {
    if (!coords) throw new Error('coords is required');

    // Feature
    if (coords.type === 'Feature' && coords.geometry !== null) return coords.geometry.coordinates;

    // Geometry
    if (coords.coordinates) return coords.coordinates;

    // Array of numbers
    if (Array.isArray(coords)) return coords;

    throw new Error('coords must be GeoJSON Feature, Geometry Object or an Array');
}

/**
 * Checks if coordinates contains a number
 *
 * @name containsNumber
 * @param {Array<any>} coordinates GeoJSON Coordinates
 * @returns {boolean} true if Array contains a number
 */
function containsNumber(coordinates) {
    if (coordinates.length > 1 && helpers.isNumber(coordinates[0]) && helpers.isNumber(coordinates[1])) {
        return true;
    }

    if (Array.isArray(coordinates[0]) && coordinates[0].length) {
        return containsNumber(coordinates[0]);
    }
    throw new Error('coordinates must only contain numbers');
}

/**
 * Enforce expectations about types of GeoJSON objects for Turf.
 *
 * @name geojsonType
 * @param {GeoJSON} value any GeoJSON object
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function geojsonType(value, type, name) {
    if (!type || !name) throw new Error('type and name required');

    if (!value || value.type !== type) {
        throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + value.type);
    }
}

/**
 * Enforce expectations about types of {@link Feature} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name featureOf
 * @param {Feature} feature a feature with an expected geometry type
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} error if value is not the expected type.
 */
function featureOf(feature, type, name) {
    if (!feature) throw new Error('No feature passed');
    if (!name) throw new Error('.featureOf() requires a name');
    if (!feature || feature.type !== 'Feature' || !feature.geometry) {
        throw new Error('Invalid input to ' + name + ', Feature with geometry required');
    }
    if (!feature.geometry || feature.geometry.type !== type) {
        throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + feature.geometry.type);
    }
}

/**
 * Enforce expectations about types of {@link FeatureCollection} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name collectionOf
 * @param {FeatureCollection} featureCollection a FeatureCollection for which features will be judged
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function collectionOf(featureCollection, type, name) {
    if (!featureCollection) throw new Error('No featureCollection passed');
    if (!name) throw new Error('.collectionOf() requires a name');
    if (!featureCollection || featureCollection.type !== 'FeatureCollection') {
        throw new Error('Invalid input to ' + name + ', FeatureCollection required');
    }
    for (var i = 0; i < featureCollection.features.length; i++) {
        var feature = featureCollection.features[i];
        if (!feature || feature.type !== 'Feature' || !feature.geometry) {
            throw new Error('Invalid input to ' + name + ', Feature with geometry required');
        }
        if (!feature.geometry || feature.geometry.type !== type) {
            throw new Error('Invalid input to ' + name + ': must be a ' + type + ', given ' + feature.geometry.type);
        }
    }
}

/**
 * Get Geometry from Feature or Geometry Object
 *
 * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
 * @returns {Geometry|null} GeoJSON Geometry Object
 * @throws {Error} if geojson is not a Feature or Geometry Object
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getGeom(point)
 * //={"type": "Point", "coordinates": [110, 40]}
 */
function getGeom(geojson) {
    if (!geojson) throw new Error('geojson is required');
    if (geojson.geometry !== undefined) return geojson.geometry;
    if (geojson.coordinates || geojson.geometries) return geojson;
    throw new Error('geojson must be a valid Feature or Geometry Object');
}

/**
 * Get Geometry Type from Feature or Geometry Object
 *
 * @throws {Error} **DEPRECATED** in v5.0.0 in favor of getType
 */
function getGeomType() {
    throw new Error('invariant.getGeomType has been deprecated in v5.0 in favor of invariant.getType');
}

/**
 * Get GeoJSON object's type, Geometry type is prioritize.
 *
 * @param {GeoJSON} geojson GeoJSON object
 * @param {string} [name="geojson"] name of the variable to display in error message
 * @returns {string} GeoJSON type
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getType(point)
 * //="Point"
 */
function getType(geojson, name) {
    if (!geojson) throw new Error((name || 'geojson') + ' is required');
    // GeoJSON Feature & GeometryCollection
    if (geojson.geometry && geojson.geometry.type) return geojson.geometry.type;
    // GeoJSON Geometry & FeatureCollection
    if (geojson.type) return geojson.type;
    throw new Error((name || 'geojson') + ' is invalid');
}

exports.getCoord = getCoord;
exports.getCoords = getCoords;
exports.containsNumber = containsNumber;
exports.geojsonType = geojsonType;
exports.featureOf = featureOf;
exports.collectionOf = collectionOf;
exports.getGeom = getGeom;
exports.getGeomType = getGeomType;
exports.getType = getType;

},{"@turf/helpers":28}],30:[function(require,module,exports){
'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var rbush = _interopDefault(require('geojson-rbush'));
var lineSegment = _interopDefault(require('@turf/line-segment'));
var invariant = require('@turf/invariant');
var meta = require('@turf/meta');
var helpers = require('@turf/helpers');

/**
 * Takes any LineString or Polygon GeoJSON and returns the intersecting point(s).
 *
 * @name lineIntersect
 * @param {Geometry|FeatureCollection|Feature<LineString|MultiLineString|Polygon|MultiPolygon>} line1 any LineString or Polygon
 * @param {Geometry|FeatureCollection|Feature<LineString|MultiLineString|Polygon|MultiPolygon>} line2 any LineString or Polygon
 * @returns {FeatureCollection<Point>} point(s) that intersect both
 * @example
 * var line1 = turf.lineString([[126, -11], [129, -21]]);
 * var line2 = turf.lineString([[123, -18], [131, -14]]);
 * var intersects = turf.lineIntersect(line1, line2);
 *
 * //addToMap
 * var addToMap = [line1, line2, intersects]
 */
function lineIntersect(line1, line2) {
    var unique = {};
    var results = [];

    // First, normalize geometries to features
    // Then, handle simple 2-vertex segments
    if (line1.type === 'LineString') line1 = helpers.feature(line1);
    if (line2.type === 'LineString') line2 = helpers.feature(line2);
    if (line1.type === 'Feature' &&
        line2.type === 'Feature' &&
        line1.geometry.type === 'LineString' &&
        line2.geometry.type === 'LineString' &&
        line1.geometry.coordinates.length === 2 &&
        line2.geometry.coordinates.length === 2) {
        var intersect = intersects(line1, line2);
        if (intersect) results.push(intersect);
        return helpers.featureCollection(results);
    }

    // Handles complex GeoJSON Geometries
    var tree = rbush();
    tree.load(lineSegment(line2));
    meta.featureEach(lineSegment(line1), function (segment) {
        meta.featureEach(tree.search(segment), function (match) {
            var intersect = intersects(segment, match);
            if (intersect) {
                // prevent duplicate points https://github.com/Turfjs/turf/issues/688
                var key = invariant.getCoords(intersect).join(',');
                if (!unique[key]) {
                    unique[key] = true;
                    results.push(intersect);
                }
            }
        });
    });
    return helpers.featureCollection(results);
}

/**
 * Find a point that intersects LineStrings with two coordinates each
 *
 * @private
 * @param {Feature<LineString>} line1 GeoJSON LineString (Must only contain 2 coordinates)
 * @param {Feature<LineString>} line2 GeoJSON LineString (Must only contain 2 coordinates)
 * @returns {Feature<Point>} intersecting GeoJSON Point
 */
function intersects(line1, line2) {
    var coords1 = invariant.getCoords(line1);
    var coords2 = invariant.getCoords(line2);
    if (coords1.length !== 2) {
        throw new Error('<intersects> line1 must only contain 2 coordinates');
    }
    if (coords2.length !== 2) {
        throw new Error('<intersects> line2 must only contain 2 coordinates');
    }
    var x1 = coords1[0][0];
    var y1 = coords1[0][1];
    var x2 = coords1[1][0];
    var y2 = coords1[1][1];
    var x3 = coords2[0][0];
    var y3 = coords2[0][1];
    var x4 = coords2[1][0];
    var y4 = coords2[1][1];
    var denom = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
    var numeA = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
    var numeB = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));

    if (denom === 0) {
        if (numeA === 0 && numeB === 0) {
            return null;
        }
        return null;
    }

    var uA = numeA / denom;
    var uB = numeB / denom;

    if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
        var x = x1 + (uA * (x2 - x1));
        var y = y1 + (uA * (y2 - y1));
        return helpers.point([x, y]);
    }
    return null;
}

module.exports = lineIntersect;
module.exports.default = lineIntersect;

},{"@turf/helpers":28,"@turf/invariant":29,"@turf/line-segment":31,"@turf/meta":34,"geojson-rbush":33}],31:[function(require,module,exports){
'use strict';

var helpers = require('@turf/helpers');
var invariant = require('@turf/invariant');
var meta = require('@turf/meta');

/**
 * Creates a {@link FeatureCollection} of 2-vertex {@link LineString} segments from a {@link LineString|(Multi)LineString} or {@link Polygon|(Multi)Polygon}.
 *
 * @name lineSegment
 * @param {Geometry|FeatureCollection|Feature<LineString|MultiLineString|MultiPolygon|Polygon>} geojson GeoJSON Polygon or LineString
 * @returns {FeatureCollection<LineString>} 2-vertex line segments
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 * var segments = turf.lineSegment(polygon);
 *
 * //addToMap
 * var addToMap = [polygon, segments]
 */
function lineSegment(geojson) {
    if (!geojson) throw new Error('geojson is required');

    var results = [];
    meta.flattenEach(geojson, function (feature) {
        lineSegmentFeature(feature, results);
    });
    return helpers.featureCollection(results);
}

/**
 * Line Segment
 *
 * @private
 * @param {Feature<LineString|Polygon>} geojson Line or polygon feature
 * @param {Array} results push to results
 * @returns {void}
 */
function lineSegmentFeature(geojson, results) {
    var coords = [];
    var geometry = geojson.geometry;
    switch (geometry.type) {
    case 'Polygon':
        coords = invariant.getCoords(geometry);
        break;
    case 'LineString':
        coords = [invariant.getCoords(geometry)];
    }
    coords.forEach(function (coord) {
        var segments = createSegments(coord, geojson.properties);
        segments.forEach(function (segment) {
            segment.id = results.length;
            results.push(segment);
        });
    });
}

/**
 * Create Segments from LineString coordinates
 *
 * @private
 * @param {LineString} coords LineString coordinates
 * @param {*} properties GeoJSON properties
 * @returns {Array<Feature<LineString>>} line segments
 */
function createSegments(coords, properties) {
    var segments = [];
    coords.reduce(function (previousCoords, currentCoords) {
        var segment = helpers.lineString([previousCoords, currentCoords], properties);
        segment.bbox = bbox(previousCoords, currentCoords);
        segments.push(segment);
        return currentCoords;
    });
    return segments;
}

/**
 * Create BBox between two coordinates (faster than @turf/bbox)
 *
 * @private
 * @param {Array<number>} coords1 Point coordinate
 * @param {Array<number>} coords2 Point coordinate
 * @returns {BBox} [west, south, east, north]
 */
function bbox(coords1, coords2) {
    var x1 = coords1[0];
    var y1 = coords1[1];
    var x2 = coords2[0];
    var y2 = coords2[1];
    var west = (x1 < x2) ? x1 : x2;
    var south = (y1 < y2) ? y1 : y2;
    var east = (x1 > x2) ? x1 : x2;
    var north = (y1 > y2) ? y1 : y2;
    return [west, south, east, north];
}

module.exports = lineSegment;
module.exports.default = lineSegment;

},{"@turf/helpers":28,"@turf/invariant":29,"@turf/meta":34}],32:[function(require,module,exports){
'use strict';

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var bearing = _interopDefault(require('@turf/bearing'));
var distance = _interopDefault(require('@turf/distance'));
var destination = _interopDefault(require('@turf/destination'));
var lineIntersects = _interopDefault(require('@turf/line-intersect'));
var meta = require('@turf/meta');
var helpers = require('@turf/helpers');
var invariant = require('@turf/invariant');

/**
 * Takes a {@link Point} and a {@link LineString} and calculates the closest Point on the (Multi)LineString.
 *
 * @name nearestPointOnLine
 * @param {Geometry|Feature<LineString|MultiLineString>} lines lines to snap to
 * @param {Geometry|Feature<Point>|number[]} pt point to snap from
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
 * @returns {Feature<Point>} closest point on the `line` to `point`. The properties object will contain three values: `index`: closest point was found on nth line part, `dist`: distance between pt and the closest point, `location`: distance along the line between start and the closest point.
 * @example
 * var line = turf.lineString([
 *     [-77.031669, 38.878605],
 *     [-77.029609, 38.881946],
 *     [-77.020339, 38.884084],
 *     [-77.025661, 38.885821],
 *     [-77.021884, 38.889563],
 *     [-77.019824, 38.892368]
 * ]);
 * var pt = turf.point([-77.037076, 38.884017]);
 *
 * var snapped = turf.nearestPointOnLine(line, pt, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [line, pt, snapped];
 * snapped.properties['marker-color'] = '#00f';
 */
function nearestPointOnLine(lines, pt, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');

    // validation
    var type = (lines.geometry) ? lines.geometry.type : lines.type;
    if (type !== 'LineString' && type !== 'MultiLineString') {
        throw new Error('lines must be LineString or MultiLineString');
    }

    var closestPt = helpers.point([Infinity, Infinity], {
        dist: Infinity
    });

    var length = 0.0;
    meta.flattenEach(lines, function (line) {
        var coords = invariant.getCoords(line);

        for (var i = 0; i < coords.length - 1; i++) {
            //start
            var start = helpers.point(coords[i]);
            start.properties.dist = distance(pt, start, options);
            //stop
            var stop = helpers.point(coords[i + 1]);
            stop.properties.dist = distance(pt, stop, options);
            // sectionLength
            var sectionLength = distance(start, stop, options);
            //perpendicular
            var heightDistance = Math.max(start.properties.dist, stop.properties.dist);
            var direction = bearing(start, stop);
            var perpendicularPt1 = destination(pt, heightDistance, direction + 90, options);
            var perpendicularPt2 = destination(pt, heightDistance, direction - 90, options);
            var intersect = lineIntersects(
                helpers.lineString([perpendicularPt1.geometry.coordinates, perpendicularPt2.geometry.coordinates]),
                helpers.lineString([start.geometry.coordinates, stop.geometry.coordinates])
            );
            var intersectPt = null;
            if (intersect.features.length > 0) {
                intersectPt = intersect.features[0];
                intersectPt.properties.dist = distance(pt, intersectPt, options);
                intersectPt.properties.location = length + distance(start, intersectPt, options);
            }

            if (start.properties.dist < closestPt.properties.dist) {
                closestPt = start;
                closestPt.properties.index = i;
                closestPt.properties.location = length;
            }
            if (stop.properties.dist < closestPt.properties.dist) {
                closestPt = stop;
                closestPt.properties.index = i + 1;
                closestPt.properties.location = length + sectionLength;
            }
            if (intersectPt && intersectPt.properties.dist < closestPt.properties.dist) {
                closestPt = intersectPt;
                closestPt.properties.index = i;
            }
            // update length
            length += sectionLength;
        }

    });

    return closestPt;
}

module.exports = nearestPointOnLine;
module.exports.default = nearestPointOnLine;

},{"@turf/bearing":25,"@turf/destination":26,"@turf/distance":27,"@turf/helpers":28,"@turf/invariant":29,"@turf/line-intersect":30,"@turf/meta":34}],33:[function(require,module,exports){
'use strict';

function quickselect(arr, k, left, right, compare) {
    quickselectStep(arr, k, left || 0, right || (arr.length - 1), compare || defaultCompare);
}

function quickselectStep(arr, k, left, right, compare) {

    while (right > left) {
        if (right - left > 600) {
            var n = right - left + 1;
            var m = k - left + 1;
            var z = Math.log(n);
            var s = 0.5 * Math.exp(2 * z / 3);
            var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            quickselectStep(arr, k, newLeft, newRight, compare);
        }

        var t = arr[k];
        var i = left;
        var j = right;

        swap(arr, left, k);
        if (compare(arr[right], t) > 0) swap(arr, left, right);

        while (i < j) {
            swap(arr, i, j);
            i++;
            j--;
            while (compare(arr[i], t) < 0) i++;
            while (compare(arr[j], t) > 0) j--;
        }

        if (compare(arr[left], t) === 0) swap(arr, left, j);
        else {
            j++;
            swap(arr, j, right);
        }

        if (j <= k) left = j + 1;
        if (k <= j) right = j - 1;
    }
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function defaultCompare(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
}

function rbush(maxEntries, format) {
    if (!(this instanceof rbush)) return new rbush(maxEntries, format);

    // max entries in a node is 9 by default; min node fill is 40% for best performance
    this._maxEntries = Math.max(4, maxEntries || 9);
    this._minEntries = Math.max(2, Math.ceil(this._maxEntries * 0.4));

    if (format) {
        this._initFormat(format);
    }

    this.clear();
}

rbush.prototype = {

    all: function () {
        return this._all(this.data, []);
    },

    search: function (bbox) {

        var node = this.data,
            result = [],
            toBBox = this.toBBox;

        if (!intersects(bbox, node)) return result;

        var nodesToSearch = [],
            i, len, child, childBBox;

        while (node) {
            for (i = 0, len = node.children.length; i < len; i++) {

                child = node.children[i];
                childBBox = node.leaf ? toBBox(child) : child;

                if (intersects(bbox, childBBox)) {
                    if (node.leaf) result.push(child);
                    else if (contains(bbox, childBBox)) this._all(child, result);
                    else nodesToSearch.push(child);
                }
            }
            node = nodesToSearch.pop();
        }

        return result;
    },

    collides: function (bbox) {

        var node = this.data,
            toBBox = this.toBBox;

        if (!intersects(bbox, node)) return false;

        var nodesToSearch = [],
            i, len, child, childBBox;

        while (node) {
            for (i = 0, len = node.children.length; i < len; i++) {

                child = node.children[i];
                childBBox = node.leaf ? toBBox(child) : child;

                if (intersects(bbox, childBBox)) {
                    if (node.leaf || contains(bbox, childBBox)) return true;
                    nodesToSearch.push(child);
                }
            }
            node = nodesToSearch.pop();
        }

        return false;
    },

    load: function (data) {
        if (!(data && data.length)) return this;

        if (data.length < this._minEntries) {
            for (var i = 0, len = data.length; i < len; i++) {
                this.insert(data[i]);
            }
            return this;
        }

        // recursively build the tree with the given data from scratch using OMT algorithm
        var node = this._build(data.slice(), 0, data.length - 1, 0);

        if (!this.data.children.length) {
            // save as is if tree is empty
            this.data = node;

        } else if (this.data.height === node.height) {
            // split root if trees have the same height
            this._splitRoot(this.data, node);

        } else {
            if (this.data.height < node.height) {
                // swap trees if inserted one is bigger
                var tmpNode = this.data;
                this.data = node;
                node = tmpNode;
            }

            // insert the small tree into the large tree at appropriate level
            this._insert(node, this.data.height - node.height - 1, true);
        }

        return this;
    },

    insert: function (item) {
        if (item) this._insert(item, this.data.height - 1);
        return this;
    },

    clear: function () {
        this.data = createNode([]);
        return this;
    },

    remove: function (item, equalsFn) {
        if (!item) return this;

        var node = this.data,
            bbox = this.toBBox(item),
            path = [],
            indexes = [],
            i, parent, index, goingUp;

        // depth-first iterative tree traversal
        while (node || path.length) {

            if (!node) { // go up
                node = path.pop();
                parent = path[path.length - 1];
                i = indexes.pop();
                goingUp = true;
            }

            if (node.leaf) { // check current node
                index = findItem(item, node.children, equalsFn);

                if (index !== -1) {
                    // item found, remove the item and condense tree upwards
                    node.children.splice(index, 1);
                    path.push(node);
                    this._condense(path);
                    return this;
                }
            }

            if (!goingUp && !node.leaf && contains(node, bbox)) { // go down
                path.push(node);
                indexes.push(i);
                i = 0;
                parent = node;
                node = node.children[0];

            } else if (parent) { // go right
                i++;
                node = parent.children[i];
                goingUp = false;

            } else node = null; // nothing found
        }

        return this;
    },

    toBBox: function (item) { return item; },

    compareMinX: compareNodeMinX,
    compareMinY: compareNodeMinY,

    toJSON: function () { return this.data; },

    fromJSON: function (data) {
        this.data = data;
        return this;
    },

    _all: function (node, result) {
        var nodesToSearch = [];
        while (node) {
            if (node.leaf) result.push.apply(result, node.children);
            else nodesToSearch.push.apply(nodesToSearch, node.children);

            node = nodesToSearch.pop();
        }
        return result;
    },

    _build: function (items, left, right, height) {

        var N = right - left + 1,
            M = this._maxEntries,
            node;

        if (N <= M) {
            // reached leaf level; return leaf
            node = createNode(items.slice(left, right + 1));
            calcBBox(node, this.toBBox);
            return node;
        }

        if (!height) {
            // target height of the bulk-loaded tree
            height = Math.ceil(Math.log(N) / Math.log(M));

            // target number of root entries to maximize storage utilization
            M = Math.ceil(N / Math.pow(M, height - 1));
        }

        node = createNode([]);
        node.leaf = false;
        node.height = height;

        // split the items into M mostly square tiles

        var N2 = Math.ceil(N / M),
            N1 = N2 * Math.ceil(Math.sqrt(M)),
            i, j, right2, right3;

        multiSelect(items, left, right, N1, this.compareMinX);

        for (i = left; i <= right; i += N1) {

            right2 = Math.min(i + N1 - 1, right);

            multiSelect(items, i, right2, N2, this.compareMinY);

            for (j = i; j <= right2; j += N2) {

                right3 = Math.min(j + N2 - 1, right2);

                // pack each entry recursively
                node.children.push(this._build(items, j, right3, height - 1));
            }
        }

        calcBBox(node, this.toBBox);

        return node;
    },

    _chooseSubtree: function (bbox, node, level, path) {

        var i, len, child, targetNode, area, enlargement, minArea, minEnlargement;

        while (true) {
            path.push(node);

            if (node.leaf || path.length - 1 === level) break;

            minArea = minEnlargement = Infinity;

            for (i = 0, len = node.children.length; i < len; i++) {
                child = node.children[i];
                area = bboxArea(child);
                enlargement = enlargedArea(bbox, child) - area;

                // choose entry with the least area enlargement
                if (enlargement < minEnlargement) {
                    minEnlargement = enlargement;
                    minArea = area < minArea ? area : minArea;
                    targetNode = child;

                } else if (enlargement === minEnlargement) {
                    // otherwise choose one with the smallest area
                    if (area < minArea) {
                        minArea = area;
                        targetNode = child;
                    }
                }
            }

            node = targetNode || node.children[0];
        }

        return node;
    },

    _insert: function (item, level, isNode) {

        var toBBox = this.toBBox,
            bbox = isNode ? item : toBBox(item),
            insertPath = [];

        // find the best node for accommodating the item, saving all nodes along the path too
        var node = this._chooseSubtree(bbox, this.data, level, insertPath);

        // put the item into the node
        node.children.push(item);
        extend(node, bbox);

        // split on node overflow; propagate upwards if necessary
        while (level >= 0) {
            if (insertPath[level].children.length > this._maxEntries) {
                this._split(insertPath, level);
                level--;
            } else break;
        }

        // adjust bboxes along the insertion path
        this._adjustParentBBoxes(bbox, insertPath, level);
    },

    // split overflowed node into two
    _split: function (insertPath, level) {

        var node = insertPath[level],
            M = node.children.length,
            m = this._minEntries;

        this._chooseSplitAxis(node, m, M);

        var splitIndex = this._chooseSplitIndex(node, m, M);

        var newNode = createNode(node.children.splice(splitIndex, node.children.length - splitIndex));
        newNode.height = node.height;
        newNode.leaf = node.leaf;

        calcBBox(node, this.toBBox);
        calcBBox(newNode, this.toBBox);

        if (level) insertPath[level - 1].children.push(newNode);
        else this._splitRoot(node, newNode);
    },

    _splitRoot: function (node, newNode) {
        // split root node
        this.data = createNode([node, newNode]);
        this.data.height = node.height + 1;
        this.data.leaf = false;
        calcBBox(this.data, this.toBBox);
    },

    _chooseSplitIndex: function (node, m, M) {

        var i, bbox1, bbox2, overlap, area, minOverlap, minArea, index;

        minOverlap = minArea = Infinity;

        for (i = m; i <= M - m; i++) {
            bbox1 = distBBox(node, 0, i, this.toBBox);
            bbox2 = distBBox(node, i, M, this.toBBox);

            overlap = intersectionArea(bbox1, bbox2);
            area = bboxArea(bbox1) + bboxArea(bbox2);

            // choose distribution with minimum overlap
            if (overlap < minOverlap) {
                minOverlap = overlap;
                index = i;

                minArea = area < minArea ? area : minArea;

            } else if (overlap === minOverlap) {
                // otherwise choose distribution with minimum area
                if (area < minArea) {
                    minArea = area;
                    index = i;
                }
            }
        }

        return index;
    },

    // sorts node children by the best axis for split
    _chooseSplitAxis: function (node, m, M) {

        var compareMinX = node.leaf ? this.compareMinX : compareNodeMinX,
            compareMinY = node.leaf ? this.compareMinY : compareNodeMinY,
            xMargin = this._allDistMargin(node, m, M, compareMinX),
            yMargin = this._allDistMargin(node, m, M, compareMinY);

        // if total distributions margin value is minimal for x, sort by minX,
        // otherwise it's already sorted by minY
        if (xMargin < yMargin) node.children.sort(compareMinX);
    },

    // total margin of all possible split distributions where each node is at least m full
    _allDistMargin: function (node, m, M, compare) {

        node.children.sort(compare);

        var toBBox = this.toBBox,
            leftBBox = distBBox(node, 0, m, toBBox),
            rightBBox = distBBox(node, M - m, M, toBBox),
            margin = bboxMargin(leftBBox) + bboxMargin(rightBBox),
            i, child;

        for (i = m; i < M - m; i++) {
            child = node.children[i];
            extend(leftBBox, node.leaf ? toBBox(child) : child);
            margin += bboxMargin(leftBBox);
        }

        for (i = M - m - 1; i >= m; i--) {
            child = node.children[i];
            extend(rightBBox, node.leaf ? toBBox(child) : child);
            margin += bboxMargin(rightBBox);
        }

        return margin;
    },

    _adjustParentBBoxes: function (bbox, path, level) {
        // adjust bboxes along the given tree path
        for (var i = level; i >= 0; i--) {
            extend(path[i], bbox);
        }
    },

    _condense: function (path) {
        // go through the path, removing empty nodes and updating bboxes
        for (var i = path.length - 1, siblings; i >= 0; i--) {
            if (path[i].children.length === 0) {
                if (i > 0) {
                    siblings = path[i - 1].children;
                    siblings.splice(siblings.indexOf(path[i]), 1);

                } else this.clear();

            } else calcBBox(path[i], this.toBBox);
        }
    },

    _initFormat: function (format) {
        // data format (minX, minY, maxX, maxY accessors)

        // uses eval-type function compilation instead of just accepting a toBBox function
        // because the algorithms are very sensitive to sorting functions performance,
        // so they should be dead simple and without inner calls

        var compareArr = ['return a', ' - b', ';'];

        this.compareMinX = new Function('a', 'b', compareArr.join(format[0]));
        this.compareMinY = new Function('a', 'b', compareArr.join(format[1]));

        this.toBBox = new Function('a',
            'return {minX: a' + format[0] +
            ', minY: a' + format[1] +
            ', maxX: a' + format[2] +
            ', maxY: a' + format[3] + '};');
    }
};

function findItem(item, items, equalsFn) {
    if (!equalsFn) return items.indexOf(item);

    for (var i = 0; i < items.length; i++) {
        if (equalsFn(item, items[i])) return i;
    }
    return -1;
}

// calculate node's bbox from bboxes of its children
function calcBBox(node, toBBox) {
    distBBox(node, 0, node.children.length, toBBox, node);
}

// min bounding rectangle of node children from k to p-1
function distBBox(node, k, p, toBBox, destNode) {
    if (!destNode) destNode = createNode(null);
    destNode.minX = Infinity;
    destNode.minY = Infinity;
    destNode.maxX = -Infinity;
    destNode.maxY = -Infinity;

    for (var i = k, child; i < p; i++) {
        child = node.children[i];
        extend(destNode, node.leaf ? toBBox(child) : child);
    }

    return destNode;
}

function extend(a, b) {
    a.minX = Math.min(a.minX, b.minX);
    a.minY = Math.min(a.minY, b.minY);
    a.maxX = Math.max(a.maxX, b.maxX);
    a.maxY = Math.max(a.maxY, b.maxY);
    return a;
}

function compareNodeMinX(a, b) { return a.minX - b.minX; }
function compareNodeMinY(a, b) { return a.minY - b.minY; }

function bboxArea(a)   { return (a.maxX - a.minX) * (a.maxY - a.minY); }
function bboxMargin(a) { return (a.maxX - a.minX) + (a.maxY - a.minY); }

function enlargedArea(a, b) {
    return (Math.max(b.maxX, a.maxX) - Math.min(b.minX, a.minX)) *
           (Math.max(b.maxY, a.maxY) - Math.min(b.minY, a.minY));
}

function intersectionArea(a, b) {
    var minX = Math.max(a.minX, b.minX),
        minY = Math.max(a.minY, b.minY),
        maxX = Math.min(a.maxX, b.maxX),
        maxY = Math.min(a.maxY, b.maxY);

    return Math.max(0, maxX - minX) *
           Math.max(0, maxY - minY);
}

function contains(a, b) {
    return a.minX <= b.minX &&
           a.minY <= b.minY &&
           b.maxX <= a.maxX &&
           b.maxY <= a.maxY;
}

function intersects(a, b) {
    return b.minX <= a.maxX &&
           b.minY <= a.maxY &&
           b.maxX >= a.minX &&
           b.maxY >= a.minY;
}

function createNode(children) {
    return {
        children: children,
        height: 1,
        leaf: true,
        minX: Infinity,
        minY: Infinity,
        maxX: -Infinity,
        maxY: -Infinity
    };
}

// sort an array so that items come in groups of n unsorted items, with groups sorted between each other;
// combines selection algorithm with binary divide & conquer approach

function multiSelect(arr, left, right, n, compare) {
    var stack = [left, right],
        mid;

    while (stack.length) {
        right = stack.pop();
        left = stack.pop();

        if (right - left <= n) continue;

        mid = left + Math.ceil((right - left) / n / 2) * n;
        quickselect(arr, mid, left, right, compare);

        stack.push(left, mid, mid, right);
    }
}

/**
 * Callback for coordEach
 *
 * @callback coordEachCallback
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * Starts at index 0.
 * @param {number} featureIndex The current index of the feature being processed.
 * @param {number} featureSubIndex The current subIndex of the feature being processed.
 */

/**
 * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
 *
 * @name coordEach
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, featureSubIndex)
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, featureSubIndex) {
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=featureSubIndex
 * });
 */
function coordEach(geojson, callback, excludeWrapCoord) {
    // Handles null Geometry -- Skips this GeoJSON
    if (geojson === null) return;
    var featureIndex, geometryIndex, j, k, l, geometry, stopG, coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        coordIndex = 0,
        isGeometryCollection,
        type = geojson.type,
        isFeatureCollection = type === 'FeatureCollection',
        isFeature = type === 'Feature',
        stop = isFeatureCollection ? geojson.features.length : 1;

    // This logic may look a little weird. The reason why it is that way
    // is because it's trying to be fast. GeoJSON supports multiple kinds
    // of objects at its root: FeatureCollection, Features, Geometries.
    // This function has the responsibility of handling all of them, and that
    // means that some of the `for` loops you see below actually just don't apply
    // to certain inputs. For instance, if you give this just a
    // Point geometry, then both loops are short-circuited and all we do
    // is gradually rename the input until it's called 'geometry'.
    //
    // This also aims to allocate as few resources as possible: just a
    // few numbers and booleans, rather than any temporary arrays as would
    // be required with the normalization approach.
    for (featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = (isFeatureCollection ? geojson.features[featureIndex].geometry :
            (isFeature ? geojson.geometry : geojson));
        isGeometryCollection = (geometryMaybeCollection) ? geometryMaybeCollection.type === 'GeometryCollection' : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (geometryIndex = 0; geometryIndex < stopG; geometryIndex++) {
            var featureSubIndex = 0;
            geometry = isGeometryCollection ?
                geometryMaybeCollection.geometries[geometryIndex] : geometryMaybeCollection;

            // Handles null Geometry -- Skips this geometry
            if (geometry === null) continue;
            coords = geometry.coordinates;
            var geomType = geometry.type;

            wrapShrink = (excludeWrapCoord && (geomType === 'Polygon' || geomType === 'MultiPolygon')) ? 1 : 0;

            switch (geomType) {
            case null:
                break;
            case 'Point':
                callback(coords, coordIndex, featureIndex, featureSubIndex);
                coordIndex++;
                featureSubIndex++;
                break;
            case 'LineString':
            case 'MultiPoint':
                for (j = 0; j < coords.length; j++) {
                    callback(coords[j], coordIndex, featureIndex, featureSubIndex);
                    coordIndex++;
                    if (geomType === 'MultiPoint') featureSubIndex++;
                }
                if (geomType === 'LineString') featureSubIndex++;
                break;
            case 'Polygon':
            case 'MultiLineString':
                for (j = 0; j < coords.length; j++) {
                    for (k = 0; k < coords[j].length - wrapShrink; k++) {
                        callback(coords[j][k], coordIndex, featureIndex, featureSubIndex);
                        coordIndex++;
                    }
                    if (geomType === 'MultiLineString') featureSubIndex++;
                }
                if (geomType === 'Polygon') featureSubIndex++;
                break;
            case 'MultiPolygon':
                for (j = 0; j < coords.length; j++) {
                    for (k = 0; k < coords[j].length; k++)
                        for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                            callback(coords[j][k][l], coordIndex, featureIndex, featureSubIndex);
                            coordIndex++;
                        }
                    featureSubIndex++;
                }
                break;
            case 'GeometryCollection':
                for (j = 0; j < geometry.geometries.length; j++)
                    coordEach(geometry.geometries[j], callback, excludeWrapCoord);
                break;
            default:
                throw new Error('Unknown Geometry Type');
            }
        }
    }
}

/**
 * Callback for coordReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback coordReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureIndex The current index of the feature being processed.
 * @param {number} featureSubIndex The current subIndex of the feature being processed.
 */

/**
 * Reduce coordinates in any GeoJSON object, similar to Array.reduce()
 *
 * @name coordReduce
 * @param {FeatureCollection|Geometry|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentCoord, coordIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordReduce(features, function (previousValue, currentCoord, coordIndex, featureIndex, featureSubIndex) {
 *   //=previousValue
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=featureSubIndex
 *   return currentCoord;
 * });
 */


/**
 * Callback for propEach
 *
 * @callback propEachCallback
 * @param {Object} currentProperties The current properties being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 */

/**
 * Iterate over properties in any GeoJSON object, similar to Array.forEach()
 *
 * @name propEach
 * @param {(FeatureCollection|Feature)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentProperties, featureIndex)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propEach(features, function (currentProperties, featureIndex) {
 *   //=currentProperties
 *   //=featureIndex
 * });
 */



/**
 * Callback for propReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback propReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {*} currentProperties The current properties being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 */

/**
 * Reduce properties in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all properties is unnecessary.
 *
 * @name propReduce
 * @param {(FeatureCollection|Feature)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentProperties, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propReduce(features, function (previousValue, currentProperties, featureIndex) {
 *   //=previousValue
 *   //=currentProperties
 *   //=featureIndex
 *   return currentProperties
 * });
 */


/**
 * Callback for featureEach
 *
 * @callback featureEachCallback
 * @param {Feature<any>} currentFeature The current feature being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 */

/**
 * Iterate over features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name featureEach
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex)
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.featureEach(features, function (currentFeature, featureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 * });
 */
function featureEach(geojson, callback) {
    if (geojson.type === 'Feature') {
        callback(geojson, 0);
    } else if (geojson.type === 'FeatureCollection') {
        for (var i = 0; i < geojson.features.length; i++) {
            callback(geojson.features[i], i);
        }
    }
}

/**
 * Callback for featureReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback featureReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name featureReduce
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.featureReduce(features, function (previousValue, currentFeature, featureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   return currentFeature
 * });
 */


/**
 * Get all coordinates from any GeoJSON object.
 *
 * @name coordAll
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @returns {Array<Array<number>>} coordinate position array
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * var coords = turf.coordAll(features);
 * //= [[26, 37], [36, 53]]
 */


/**
 * Callback for geomEach
 *
 * @callback geomEachCallback
 * @param {Geometry} currentGeometry The current geometry being processed.
 * @param {number} currentIndex The index of the current element being processed in the
 * array. Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} currentProperties The current feature properties being processed.
 */

/**
 * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
 *
 * @name geomEach
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentGeometry, featureIndex, currentProperties)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomEach(features, function (currentGeometry, featureIndex, currentProperties) {
 *   //=currentGeometry
 *   //=featureIndex
 *   //=currentProperties
 * });
 */


/**
 * Callback for geomReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback geomReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Geometry} currentGeometry The current Feature being processed.
 * @param {number} currentIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {Object} currentProperties The current feature properties being processed.
 */

/**
 * Reduce geometry in any GeoJSON object, similar to Array.reduce().
 *
 * @name geomReduce
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentGeometry, featureIndex, currentProperties)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomReduce(features, function (previousValue, currentGeometry, featureIndex, currentProperties) {
 *   //=previousValue
 *   //=currentGeometry
 *   //=featureIndex
 *   //=currentProperties
 *   return currentGeometry
 * });
 */


/**
 * Callback for flattenEach
 *
 * @callback flattenEachCallback
 * @param {Feature} currentFeature The current flattened feature being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array. Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureSubIndex The subindex of the current element being processed in the
 * array. Starts at index 0 and increases if the flattened feature was a multi-geometry.
 */

/**
 * Iterate over flattened features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name flattenEach
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex, featureSubIndex)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenEach(features, function (currentFeature, featureIndex, featureSubIndex) {
 *   //=currentFeature
 *   //=featureIndex
 *   //=featureSubIndex
 * });
 */


/**
 * Callback for flattenReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback flattenReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The index of the current element being processed in the
 * array.Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureSubIndex The subindex of the current element being processed in the
 * array. Starts at index 0 and increases if the flattened feature was a multi-geometry.
 */

/**
 * Reduce flattened features in any GeoJSON object, similar to Array.reduce().
 *
 * @name flattenReduce
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex, featureSubIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenReduce(features, function (previousValue, currentFeature, featureIndex, featureSubIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   //=featureSubIndex
 *   return currentFeature
 * });
 */


/**
 * Callback for segmentEach
 *
 * @callback segmentEachCallback
 * @param {Feature<LineString>} currentSegment The current segment being processed.
 * @param {number} featureIndex The featureIndex currently being processed, starts at index 0.
 * @param {number} featureSubIndex The featureSubIndex currently being processed, starts at index 0.
 * @param {number} segmentIndex The segmentIndex currently being processed, starts at index 0.
 * @returns {void}
 */

/**
 * Iterate over 2-vertex line segment in any GeoJSON object, similar to Array.forEach()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON
 * @param {Function} callback a method that takes (currentSegment, featureIndex, featureSubIndex)
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentEach(polygon, function (currentSegment, featureIndex, featureSubIndex, segmentIndex) {
 *   //= currentSegment
 *   //= featureIndex
 *   //= featureSubIndex
 *   //= segmentIndex
 * });
 *
 * // Calculate the total number of segments
 * var total = 0;
 * turf.segmentEach(polygon, function () {
 *     total++;
 * });
 */


/**
 * Callback for segmentReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback segmentReduceCallback
 * @param {*} [previousValue] The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} [currentSegment] The current segment being processed.
 * @param {number} featureIndex The featureIndex currently being processed, starts at index 0.
 * @param {number} featureSubIndex The featureSubIndex currently being processed, starts at index 0.
 * @param {number} segmentIndex The segmentIndex currently being processed, starts at index 0.
 */

/**
 * Reduce 2-vertex line segment in any GeoJSON object, similar to Array.reduce()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {(FeatureCollection|Feature|Geometry)} geojson any GeoJSON
 * @param {Function} callback a method that takes (previousValue, currentSegment, currentIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentReduce(polygon, function (previousSegment, currentSegment, featureIndex, featureSubIndex, segmentIndex) {
 *   //= previousSegment
 *   //= currentSegment
 *   //= featureIndex
 *   //= featureSubIndex
 *   //= segmentInex
 *   return currentSegment
 * });
 *
 * // Calculate the total number of segments
 * var initialValue = 0
 * var total = turf.segmentReduce(polygon, function (previousValue) {
 *     previousValue++;
 *     return previousValue;
 * }, initialValue);
 */


/**
 * Callback for lineEach
 *
 * @callback lineEachCallback
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
 * @param {number} lineIndex The index of the current element being processed in the array, starts at index 0.
 * @param {number} lineSubIndex The sub-index of the current line being processed at index 0
 */

/**
 * Iterate over line or ring coordinates in LineString, Polygon, MultiLineString, MultiPolygon Features or Geometries,
 * similar to Array.forEach.
 *
 * @name lineEach
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (currentLine, lineIndex, lineSubIndex)
 * @example
 * var mtLn = turf.multiLineString([
 *   turf.lineString([[26, 37], [35, 45]]),
 *   turf.lineString([[36, 53], [38, 50], [41, 55]])
 * ]);
 *
 * turf.lineEach(mtLn, function (currentLine, lineIndex) {
 *   //=currentLine
 *   //=lineIndex
 * });
 */


/**
 * Callback for lineReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback lineReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
 * @param {number} lineIndex The index of the current element being processed in the
 * array. Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} lineSubIndex The sub-index of the current line being processed at index 0
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name lineReduce
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var mtp = turf.multiPolygon([
 *   turf.polygon([[[12,48],[2,41],[24,38],[12,48]], [[9,44],[13,41],[13,45],[9,44]]]),
 *   turf.polygon([[[5, 5], [0, 0], [2, 2], [4, 4], [5, 5]]])
 * ]);
 *
 * turf.lineReduce(mtp, function (previousValue, currentLine, lineIndex, lineSubIndex) {
 *   //=previousValue
 *   //=currentLine
 *   //=lineIndex
 *   //=lineSubIndex
 *   return currentLine
 * }, 2);
 */

/**
 * GeoJSON implementation of [RBush](https://github.com/mourner/rbush#rbush) spatial index.
 *
 * @name rbush
 * @param {number} [maxEntries=9] defines the maximum number of entries in a tree node. 9 (used by default) is a
 * reasonable choice for most applications. Higher value means faster insertion and slower search, and vice versa.
 * @returns {RBush} GeoJSON RBush
 * @example
 * import geojsonRbush from 'geojson-rbush';
 * var tree = geojsonRbush();
 */
function geojsonRbush(maxEntries) {
    var tree = rbush(maxEntries);
    /**
     * [insert](https://github.com/mourner/rbush#data-format)
     *
     * @param {Feature<any>} feature insert single GeoJSON Feature
     * @returns {RBush} GeoJSON RBush
     * @example
     * var polygon = {
     *   "type": "Feature",
     *   "properties": {},
     *   "geometry": {
     *     "type": "Polygon",
     *     "coordinates": [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]
     *   }
     * }
     * tree.insert(polygon)
     */
    tree.insert = function (feature) {
        if (Array.isArray(feature)) {
            var bbox = feature;
            feature = bboxPolygon(bbox);
            feature.bbox = bbox;
        } else {
            feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
        }
        return rbush.prototype.insert.call(this, feature);
    };

    /**
     * [load](https://github.com/mourner/rbush#bulk-inserting-data)
     *
     * @param {BBox[]|FeatureCollection<any>} features load entire GeoJSON FeatureCollection
     * @returns {RBush} GeoJSON RBush
     * @example
     * var polygons = {
     *   "type": "FeatureCollection",
     *   "features": [
     *     {
     *       "type": "Feature",
     *       "properties": {},
     *       "geometry": {
     *         "type": "Polygon",
     *         "coordinates": [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]
     *       }
     *     },
     *     {
     *       "type": "Feature",
     *       "properties": {},
     *       "geometry": {
     *         "type": "Polygon",
     *         "coordinates": [[[-93, 32], [-83, 32], [-83, 39], [-93, 39], [-93, 32]]]
     *       }
     *     }
     *   ]
     * }
     * tree.load(polygons)
     */
    tree.load = function (features) {
        var load = [];
        // Load an Array of BBox
        if (Array.isArray(features)) {
            features.forEach(function (bbox) {
                var feature = bboxPolygon(bbox);
                feature.bbox = bbox;
                load.push(feature);
            });
        } else {
            // Load FeatureCollection
            featureEach(features, function (feature) {
                feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
                load.push(feature);
            });
        }
        return rbush.prototype.load.call(this, load);
    };

    /**
     * [remove](https://github.com/mourner/rbush#removing-data)
     *
     * @param {BBox|Feature<any>} feature remove single GeoJSON Feature
     * @returns {RBush} GeoJSON RBush
     * @example
     * var polygon = {
     *   "type": "Feature",
     *   "properties": {},
     *   "geometry": {
     *     "type": "Polygon",
     *     "coordinates": [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]
     *   }
     * }
     * tree.remove(polygon)
     */
    tree.remove = function (feature) {
        if (Array.isArray(feature)) {
            var bbox = feature;
            feature = bboxPolygon(bbox);
            feature.bbox = bbox;
        }
        return rbush.prototype.remove.call(this, feature);
    };

    /**
     * [clear](https://github.com/mourner/rbush#removing-data)
     *
     * @returns {RBush} GeoJSON Rbush
     * @example
     * tree.clear()
     */
    tree.clear = function () {
        return rbush.prototype.clear.call(this);
    };

    /**
     * [search](https://github.com/mourner/rbush#search)
     *
     * @param {BBox|FeatureCollection|Feature<any>} geojson search with GeoJSON
     * @returns {FeatureCollection<any>} all features that intersects with the given GeoJSON.
     * @example
     * var polygon = {
     *   "type": "Feature",
     *   "properties": {},
     *   "geometry": {
     *     "type": "Polygon",
     *     "coordinates": [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]
     *   }
     * }
     * tree.search(polygon)
     */
    tree.search = function (geojson) {
        var features = rbush.prototype.search.call(this, this.toBBox(geojson));
        return {
            type: 'FeatureCollection',
            features: features
        };
    };

    /**
     * [collides](https://github.com/mourner/rbush#collisions)
     *
     * @param {BBox|FeatureCollection|Feature<any>} geojson collides with GeoJSON
     * @returns {boolean} true if there are any items intersecting the given GeoJSON, otherwise false.
     * @example
     * var polygon = {
     *   "type": "Feature",
     *   "properties": {},
     *   "geometry": {
     *     "type": "Polygon",
     *     "coordinates": [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]
     *   }
     * }
     * tree.collides(polygon)
     */
    tree.collides = function (geojson) {
        return rbush.prototype.collides.call(this, this.toBBox(geojson));
    };

    /**
     * [all](https://github.com/mourner/rbush#search)
     *
     * @returns {FeatureCollection<any>} all the features in RBush
     * @example
     * tree.all()
     * //=FeatureCollection
     */
    tree.all = function () {
        var features = rbush.prototype.all.call(this);
        return {
            type: 'FeatureCollection',
            features: features
        };
    };

    /**
     * [toJSON](https://github.com/mourner/rbush#export-and-import)
     *
     * @returns {any} export data as JSON object
     * @example
     * var exported = tree.toJSON()
     * //=JSON object
     */
    tree.toJSON = function () {
        return rbush.prototype.toJSON.call(this);
    };

    /**
     * [fromJSON](https://github.com/mourner/rbush#export-and-import)
     *
     * @param {any} json import previously exported data
     * @returns {RBush} GeoJSON RBush
     * @example
     * var exported = {
     *   "children": [
     *     {
     *       "type": "Feature",
     *       "geometry": {
     *         "type": "Point",
     *         "coordinates": [110, 50]
     *       },
     *       "properties": {},
     *       "bbox": [110, 50, 110, 50]
     *     }
     *   ],
     *   "height": 1,
     *   "leaf": true,
     *   "minX": 110,
     *   "minY": 50,
     *   "maxX": 110,
     *   "maxY": 50
     * }
     * tree.fromJSON(exported)
     */
    tree.fromJSON = function (json) {
        return rbush.prototype.fromJSON.call(this, json);
    };

    /**
     * Converts GeoJSON to {minX, minY, maxX, maxY} schema
     *
     * @private
     * @param {BBox|FeatureCollectio|Feature<any>} geojson feature(s) to retrieve BBox from
     * @returns {Object} converted to {minX, minY, maxX, maxY}
     */
    tree.toBBox = function (geojson) {
        var bbox;
        if (geojson.bbox) bbox = geojson.bbox;
        else if (Array.isArray(geojson) && geojson.length === 4) bbox = geojson;
        else bbox = turfBBox(geojson);

        return {
            minX: bbox[0],
            minY: bbox[1],
            maxX: bbox[2],
            maxY: bbox[3]
        };
    };
    return tree;
}

/**
 * Takes a bbox and returns an equivalent {@link Polygon|polygon}.
 *
 * @private
 * @name bboxPolygon
 * @param {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
 * @returns {Feature<Polygon>} a Polygon representation of the bounding box
 * @example
 * var bbox = [0, 0, 10, 10];
 *
 * var poly = turf.bboxPolygon(bbox);
 *
 * //addToMap
 * var addToMap = [poly]
 */
function bboxPolygon(bbox) {
    var lowLeft = [bbox[0], bbox[1]];
    var topLeft = [bbox[0], bbox[3]];
    var topRight = [bbox[2], bbox[3]];
    var lowRight = [bbox[2], bbox[1]];
    var coordinates = [[lowLeft, lowRight, topRight, topLeft, lowLeft]];

    return {
        type: 'Feature',
        bbox: bbox,
        properties: {},
        geometry: {
            type: 'Polygon',
            coordinates: coordinates
        }
    };
}

/**
 * Takes a set of features, calculates the bbox of all input features, and returns a bounding box.
 *
 * @private
 * @name bbox
 * @param {FeatureCollection|Feature<any>} geojson input features
 * @returns {Array<number>} bbox extent in [minX, minY, maxX, maxY] order
 * @example
 * var line = turf.lineString([[-74, 40], [-78, 42], [-82, 35]]);
 * var bbox = turf.bbox(line);
 * var bboxPolygon = turf.bboxPolygon(bbox);
 *
 * //addToMap
 * var addToMap = [line, bboxPolygon]
 */
function turfBBox(geojson) {
    var bbox = [Infinity, Infinity, -Infinity, -Infinity];
    coordEach(geojson, function (coord) {
        if (bbox[0] > coord[0]) bbox[0] = coord[0];
        if (bbox[1] > coord[1]) bbox[1] = coord[1];
        if (bbox[2] < coord[0]) bbox[2] = coord[0];
        if (bbox[3] < coord[1]) bbox[3] = coord[1];
    });
    return bbox;
}

module.exports = geojsonRbush;
module.exports.default = geojsonRbush;

},{}],34:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

var helpers = require('@turf/helpers');

/**
 * Callback for coordEach
 *
 * @callback coordEachCallback
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
 *
 * @name coordEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function coordEach(geojson, callback, excludeWrapCoord) {
    // Handles null Geometry -- Skips this GeoJSON
    if (geojson === null) return;
    var j, k, l, geometry, stopG, coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        coordIndex = 0,
        isGeometryCollection,
        type = geojson.type,
        isFeatureCollection = type === 'FeatureCollection',
        isFeature = type === 'Feature',
        stop = isFeatureCollection ? geojson.features.length : 1;

    // This logic may look a little weird. The reason why it is that way
    // is because it's trying to be fast. GeoJSON supports multiple kinds
    // of objects at its root: FeatureCollection, Features, Geometries.
    // This function has the responsibility of handling all of them, and that
    // means that some of the `for` loops you see below actually just don't apply
    // to certain inputs. For instance, if you give this just a
    // Point geometry, then both loops are short-circuited and all we do
    // is gradually rename the input until it's called 'geometry'.
    //
    // This also aims to allocate as few resources as possible: just a
    // few numbers and booleans, rather than any temporary arrays as would
    // be required with the normalization approach.
    for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = (isFeatureCollection ? geojson.features[featureIndex].geometry :
            (isFeature ? geojson.geometry : geojson));
        isGeometryCollection = (geometryMaybeCollection) ? geometryMaybeCollection.type === 'GeometryCollection' : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
            var multiFeatureIndex = 0;
            var geometryIndex = 0;
            geometry = isGeometryCollection ?
                geometryMaybeCollection.geometries[geomIndex] : geometryMaybeCollection;

            // Handles null Geometry -- Skips this geometry
            if (geometry === null) continue;
            coords = geometry.coordinates;
            var geomType = geometry.type;

            wrapShrink = (excludeWrapCoord && (geomType === 'Polygon' || geomType === 'MultiPolygon')) ? 1 : 0;

            switch (geomType) {
            case null:
                break;
            case 'Point':
                if (callback(coords, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                coordIndex++;
                multiFeatureIndex++;
                break;
            case 'LineString':
            case 'MultiPoint':
                for (j = 0; j < coords.length; j++) {
                    if (callback(coords[j], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                    coordIndex++;
                    if (geomType === 'MultiPoint') multiFeatureIndex++;
                }
                if (geomType === 'LineString') multiFeatureIndex++;
                break;
            case 'Polygon':
            case 'MultiLineString':
                for (j = 0; j < coords.length; j++) {
                    for (k = 0; k < coords[j].length - wrapShrink; k++) {
                        if (callback(coords[j][k], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                        coordIndex++;
                    }
                    if (geomType === 'MultiLineString') multiFeatureIndex++;
                    if (geomType === 'Polygon') geometryIndex++;
                }
                if (geomType === 'Polygon') multiFeatureIndex++;
                break;
            case 'MultiPolygon':
                for (j = 0; j < coords.length; j++) {
                    if (geomType === 'MultiPolygon') geometryIndex = 0;
                    for (k = 0; k < coords[j].length; k++) {
                        for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                            if (callback(coords[j][k][l], coordIndex, featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
                            coordIndex++;
                        }
                        geometryIndex++;
                    }
                    multiFeatureIndex++;
                }
                break;
            case 'GeometryCollection':
                for (j = 0; j < geometry.geometries.length; j++)
                    if (coordEach(geometry.geometries[j], callback, excludeWrapCoord) === false) return false;
                break;
            default:
                throw new Error('Unknown Geometry Type');
            }
        }
    }
}

/**
 * Callback for coordReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback coordReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Array<number>} currentCoord The current coordinate being processed.
 * @param {number} coordIndex The current index of the coordinate being processed.
 * Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 */

/**
 * Reduce coordinates in any GeoJSON object, similar to Array.reduce()
 *
 * @name coordReduce
 * @param {FeatureCollection|Geometry|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentCoord, coordIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.coordReduce(features, function (previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentCoord
 *   //=coordIndex
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentCoord;
 * });
 */
function coordReduce(geojson, callback, initialValue, excludeWrapCoord) {
    var previousValue = initialValue;
    coordEach(geojson, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
        if (coordIndex === 0 && initialValue === undefined) previousValue = currentCoord;
        else previousValue = callback(previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex);
    }, excludeWrapCoord);
    return previousValue;
}

/**
 * Callback for propEach
 *
 * @callback propEachCallback
 * @param {Object} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over properties in any GeoJSON object, similar to Array.forEach()
 *
 * @name propEach
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentProperties, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propEach(features, function (currentProperties, featureIndex) {
 *   //=currentProperties
 *   //=featureIndex
 * });
 */
function propEach(geojson, callback) {
    var i;
    switch (geojson.type) {
    case 'FeatureCollection':
        for (i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i].properties, i) === false) break;
        }
        break;
    case 'Feature':
        callback(geojson.properties, 0);
        break;
    }
}


/**
 * Callback for propReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback propReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {*} currentProperties The current Properties being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce properties in any GeoJSON object into a single value,
 * similar to how Array.reduce works. However, in this case we lazily run
 * the reduction, so an array of all properties is unnecessary.
 *
 * @name propReduce
 * @param {FeatureCollection|Feature} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentProperties, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.propReduce(features, function (previousValue, currentProperties, featureIndex) {
 *   //=previousValue
 *   //=currentProperties
 *   //=featureIndex
 *   return currentProperties
 * });
 */
function propReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    propEach(geojson, function (currentProperties, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentProperties;
        else previousValue = callback(previousValue, currentProperties, featureIndex);
    });
    return previousValue;
}

/**
 * Callback for featureEach
 *
 * @callback featureEachCallback
 * @param {Feature<any>} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Iterate over features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name featureEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.featureEach(features, function (currentFeature, featureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 * });
 */
function featureEach(geojson, callback) {
    if (geojson.type === 'Feature') {
        callback(geojson, 0);
    } else if (geojson.type === 'FeatureCollection') {
        for (var i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i], i) === false) break;
        }
    }
}

/**
 * Callback for featureReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback featureReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name featureReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {"foo": "bar"}),
 *   turf.point([36, 53], {"hello": "world"})
 * ]);
 *
 * turf.featureReduce(features, function (previousValue, currentFeature, featureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   return currentFeature
 * });
 */
function featureReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    featureEach(geojson, function (currentFeature, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentFeature;
        else previousValue = callback(previousValue, currentFeature, featureIndex);
    });
    return previousValue;
}

/**
 * Get all coordinates from any GeoJSON object.
 *
 * @name coordAll
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @returns {Array<Array<number>>} coordinate position array
 * @example
 * var features = turf.featureCollection([
 *   turf.point([26, 37], {foo: 'bar'}),
 *   turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * var coords = turf.coordAll(features);
 * //= [[26, 37], [36, 53]]
 */
function coordAll(geojson) {
    var coords = [];
    coordEach(geojson, function (coord) {
        coords.push(coord);
    });
    return coords;
}

/**
 * Callback for geomEach
 *
 * @callback geomEachCallback
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
 *
 * @name geomEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @returns {void}
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomEach(features, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 * });
 */
function geomEach(geojson, callback) {
    var i, j, g, geometry, stopG,
        geometryMaybeCollection,
        isGeometryCollection,
        featureProperties,
        featureBBox,
        featureId,
        featureIndex = 0,
        isFeatureCollection = geojson.type === 'FeatureCollection',
        isFeature = geojson.type === 'Feature',
        stop = isFeatureCollection ? geojson.features.length : 1;

    // This logic may look a little weird. The reason why it is that way
    // is because it's trying to be fast. GeoJSON supports multiple kinds
    // of objects at its root: FeatureCollection, Features, Geometries.
    // This function has the responsibility of handling all of them, and that
    // means that some of the `for` loops you see below actually just don't apply
    // to certain inputs. For instance, if you give this just a
    // Point geometry, then both loops are short-circuited and all we do
    // is gradually rename the input until it's called 'geometry'.
    //
    // This also aims to allocate as few resources as possible: just a
    // few numbers and booleans, rather than any temporary arrays as would
    // be required with the normalization approach.
    for (i = 0; i < stop; i++) {

        geometryMaybeCollection = (isFeatureCollection ? geojson.features[i].geometry :
            (isFeature ? geojson.geometry : geojson));
        featureProperties = (isFeatureCollection ? geojson.features[i].properties :
            (isFeature ? geojson.properties : {}));
        featureBBox = (isFeatureCollection ? geojson.features[i].bbox :
            (isFeature ? geojson.bbox : undefined));
        featureId = (isFeatureCollection ? geojson.features[i].id :
            (isFeature ? geojson.id : undefined));
        isGeometryCollection = (geometryMaybeCollection) ? geometryMaybeCollection.type === 'GeometryCollection' : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;

        for (g = 0; g < stopG; g++) {
            geometry = isGeometryCollection ?
                geometryMaybeCollection.geometries[g] : geometryMaybeCollection;

            // Handle null Geometry
            if (geometry === null) {
                if (callback(null, featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                continue;
            }
            switch (geometry.type) {
            case 'Point':
            case 'LineString':
            case 'MultiPoint':
            case 'Polygon':
            case 'MultiLineString':
            case 'MultiPolygon': {
                if (callback(geometry, featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                break;
            }
            case 'GeometryCollection': {
                for (j = 0; j < geometry.geometries.length; j++) {
                    if (callback(geometry.geometries[j], featureIndex, featureProperties, featureBBox, featureId) === false) return false;
                }
                break;
            }
            default:
                throw new Error('Unknown Geometry Type');
            }
        }
        // Only increase `featureIndex` per each feature
        featureIndex++;
    }
}

/**
 * Callback for geomReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback geomReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Geometry} currentGeometry The current Geometry being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {Object} featureProperties The current Feature Properties being processed.
 * @param {Array<number>} featureBBox The current Feature BBox being processed.
 * @param {number|string} featureId The current Feature Id being processed.
 */

/**
 * Reduce geometry in any GeoJSON object, similar to Array.reduce().
 *
 * @name geomReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.point([36, 53], {hello: 'world'})
 * ]);
 *
 * turf.geomReduce(features, function (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
 *   //=previousValue
 *   //=currentGeometry
 *   //=featureIndex
 *   //=featureProperties
 *   //=featureBBox
 *   //=featureId
 *   return currentGeometry
 * });
 */
function geomReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    geomEach(geojson, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentGeometry;
        else previousValue = callback(previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId);
    });
    return previousValue;
}

/**
 * Callback for flattenEach
 *
 * @callback flattenEachCallback
 * @param {Feature} currentFeature The current flattened feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Iterate over flattened features in any GeoJSON object, similar to
 * Array.forEach.
 *
 * @name flattenEach
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (currentFeature, featureIndex, multiFeatureIndex)
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenEach(features, function (currentFeature, featureIndex, multiFeatureIndex) {
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 * });
 */
function flattenEach(geojson, callback) {
    geomEach(geojson, function (geometry, featureIndex, properties, bbox, id) {
        // Callback for single geometry
        var type = (geometry === null) ? null : geometry.type;
        switch (type) {
        case null:
        case 'Point':
        case 'LineString':
        case 'Polygon':
            if (callback(helpers.feature(geometry, properties, {bbox: bbox, id: id}), featureIndex, 0) === false) return false;
            return;
        }

        var geomType;

        // Callback for multi-geometry
        switch (type) {
        case 'MultiPoint':
            geomType = 'Point';
            break;
        case 'MultiLineString':
            geomType = 'LineString';
            break;
        case 'MultiPolygon':
            geomType = 'Polygon';
            break;
        }

        for (var multiFeatureIndex = 0; multiFeatureIndex < geometry.coordinates.length; multiFeatureIndex++) {
            var coordinate = geometry.coordinates[multiFeatureIndex];
            var geom = {
                type: geomType,
                coordinates: coordinate
            };
            if (callback(helpers.feature(geom, properties), featureIndex, multiFeatureIndex) === false) return false;
        }
    });
}

/**
 * Callback for flattenReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback flattenReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature} currentFeature The current Feature being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 */

/**
 * Reduce flattened features in any GeoJSON object, similar to Array.reduce().
 *
 * @name flattenReduce
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
 * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex, multiFeatureIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var features = turf.featureCollection([
 *     turf.point([26, 37], {foo: 'bar'}),
 *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
 * ]);
 *
 * turf.flattenReduce(features, function (previousValue, currentFeature, featureIndex, multiFeatureIndex) {
 *   //=previousValue
 *   //=currentFeature
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   return currentFeature
 * });
 */
function flattenReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    flattenEach(geojson, function (currentFeature, featureIndex, multiFeatureIndex) {
        if (featureIndex === 0 && multiFeatureIndex === 0 && initialValue === undefined) previousValue = currentFeature;
        else previousValue = callback(previousValue, currentFeature, featureIndex, multiFeatureIndex);
    });
    return previousValue;
}

/**
 * Callback for segmentEach
 *
 * @callback segmentEachCallback
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 * @returns {void}
 */

/**
 * Iterate over 2-vertex line segment in any GeoJSON object, similar to Array.forEach()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex)
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentEach(polygon, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //=currentSegment
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   //=segmentIndex
 * });
 *
 * // Calculate the total number of segments
 * var total = 0;
 * turf.segmentEach(polygon, function () {
 *     total++;
 * });
 */
function segmentEach(geojson, callback) {
    flattenEach(geojson, function (feature$$1, featureIndex, multiFeatureIndex) {
        var segmentIndex = 0;

        // Exclude null Geometries
        if (!feature$$1.geometry) return;
        // (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
        var type = feature$$1.geometry.type;
        if (type === 'Point' || type === 'MultiPoint') return;

        // Generate 2-vertex line segments
        var previousCoords;
        if (coordEach(feature$$1, function (currentCoord, coordIndex, featureIndexCoord, mutliPartIndexCoord, geometryIndex) {
            // Simulating a meta.coordReduce() since `reduce` operations cannot be stopped by returning `false`
            if (previousCoords === undefined) {
                previousCoords = currentCoord;
                return;
            }
            var currentSegment = helpers.lineString([previousCoords, currentCoord], feature$$1.properties);
            if (callback(currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) === false) return false;
            segmentIndex++;
            previousCoords = currentCoord;
        }) === false) return false;
    });
}

/**
 * Callback for segmentReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback segmentReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentSegment The current Segment being processed.
 * @param {number} featureIndex The current index of the Feature being processed.
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
 * @param {number} geometryIndex The current index of the Geometry being processed.
 * @param {number} segmentIndex The current index of the Segment being processed.
 */

/**
 * Reduce 2-vertex line segment in any GeoJSON object, similar to Array.reduce()
 * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
 * @param {Function} callback a method that takes (previousValue, currentSegment, currentIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {void}
 * @example
 * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
 *
 * // Iterate over GeoJSON by 2-vertex segments
 * turf.segmentReduce(polygon, function (previousSegment, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
 *   //= previousSegment
 *   //= currentSegment
 *   //= featureIndex
 *   //= multiFeatureIndex
 *   //= geometryIndex
 *   //= segmentInex
 *   return currentSegment
 * });
 *
 * // Calculate the total number of segments
 * var initialValue = 0
 * var total = turf.segmentReduce(polygon, function (previousValue) {
 *     previousValue++;
 *     return previousValue;
 * }, initialValue);
 */
function segmentReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    var started = false;
    segmentEach(geojson, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
        if (started === false && initialValue === undefined) previousValue = currentSegment;
        else previousValue = callback(previousValue, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex);
        started = true;
    });
    return previousValue;
}

/**
 * Callback for lineEach
 *
 * @callback lineEachCallback
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Iterate over line or ring coordinates in LineString, Polygon, MultiLineString, MultiPolygon Features or Geometries,
 * similar to Array.forEach.
 *
 * @name lineEach
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @example
 * var multiLine = turf.multiLineString([
 *   [[26, 37], [35, 45]],
 *   [[36, 53], [38, 50], [41, 55]]
 * ]);
 *
 * turf.lineEach(multiLine, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 * });
 */
function lineEach(geojson, callback) {
    // validation
    if (!geojson) throw new Error('geojson is required');

    flattenEach(geojson, function (feature$$1, featureIndex, multiFeatureIndex) {
        if (feature$$1.geometry === null) return;
        var type = feature$$1.geometry.type;
        var coords = feature$$1.geometry.coordinates;
        switch (type) {
        case 'LineString':
            if (callback(feature$$1, featureIndex, multiFeatureIndex, 0, 0) === false) return false;
            break;
        case 'Polygon':
            for (var geometryIndex = 0; geometryIndex < coords.length; geometryIndex++) {
                if (callback(helpers.lineString(coords[geometryIndex], feature$$1.properties), featureIndex, multiFeatureIndex, geometryIndex) === false) return false;
            }
            break;
        }
    });
}

/**
 * Callback for lineReduce
 *
 * The first time the callback function is called, the values provided as arguments depend
 * on whether the reduce method has an initialValue argument.
 *
 * If an initialValue is provided to the reduce method:
 *  - The previousValue argument is initialValue.
 *  - The currentValue argument is the value of the first element present in the array.
 *
 * If an initialValue is not provided:
 *  - The previousValue argument is the value of the first element present in the array.
 *  - The currentValue argument is the value of the second element present in the array.
 *
 * @callback lineReduceCallback
 * @param {*} previousValue The accumulated value previously returned in the last invocation
 * of the callback, or initialValue, if supplied.
 * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
 * @param {number} featureIndex The current index of the Feature being processed
 * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
 * @param {number} geometryIndex The current index of the Geometry being processed
 */

/**
 * Reduce features in any GeoJSON object, similar to Array.reduce().
 *
 * @name lineReduce
 * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
 * @param {Function} callback a method that takes (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex)
 * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
 * @returns {*} The value that results from the reduction.
 * @example
 * var multiPoly = turf.multiPolygon([
 *   turf.polygon([[[12,48],[2,41],[24,38],[12,48]], [[9,44],[13,41],[13,45],[9,44]]]),
 *   turf.polygon([[[5, 5], [0, 0], [2, 2], [4, 4], [5, 5]]])
 * ]);
 *
 * turf.lineReduce(multiPoly, function (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
 *   //=previousValue
 *   //=currentLine
 *   //=featureIndex
 *   //=multiFeatureIndex
 *   //=geometryIndex
 *   return currentLine
 * });
 */
function lineReduce(geojson, callback, initialValue) {
    var previousValue = initialValue;
    lineEach(geojson, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
        if (featureIndex === 0 && initialValue === undefined) previousValue = currentLine;
        else previousValue = callback(previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex);
    });
    return previousValue;
}

/**
 * Finds a particular 2-vertex LineString Segment from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 * Point & MultiPoint will always return null.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.segmentIndex=0] Segment Index
 * @param {Object} [options.properties={}] Translate Properties to output LineString
 * @param {BBox} [options.bbox={}] Translate BBox to output LineString
 * @param {number|string} [options.id={}] Translate Id to output LineString
 * @returns {Feature<LineString>} 2-vertex GeoJSON Feature LineString
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findSegment(multiLine);
 * // => Feature<LineString<[[10, 10], [50, 30]]>>
 *
 * // First Segment of 2nd Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: 1});
 * // => Feature<LineString<[[-10, -10], [-50, -30]]>>
 *
 * // Last Segment of Last Multi Feature
 * turf.findSegment(multiLine, {multiFeatureIndex: -1, segmentIndex: -1});
 * // => Feature<LineString<[[-50, -30], [-30, -40]]>>
 */
function findSegment(geojson, options) {
    // Optional Parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var featureIndex = options.featureIndex || 0;
    var multiFeatureIndex = options.multiFeatureIndex || 0;
    var geometryIndex = options.geometryIndex || 0;
    var segmentIndex = options.segmentIndex || 0;

    // Find FeatureIndex
    var properties = options.properties;
    var geometry;

    switch (geojson.type) {
    case 'FeatureCollection':
        if (featureIndex < 0) featureIndex = geojson.features.length + featureIndex;
        properties = properties || geojson.features[featureIndex].properties;
        geometry = geojson.features[featureIndex].geometry;
        break;
    case 'Feature':
        properties = properties || geojson.properties;
        geometry = geojson.geometry;
        break;
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
    case 'Polygon':
    case 'MultiLineString':
    case 'MultiPolygon':
        geometry = geojson;
        break;
    default:
        throw new Error('geojson is invalid');
    }

    // Find SegmentIndex
    if (geometry === null) return null;
    var coords = geometry.coordinates;
    switch (geometry.type) {
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
        if (segmentIndex < 0) segmentIndex = coords.length + segmentIndex - 1;
        return helpers.lineString([coords[segmentIndex], coords[segmentIndex + 1]], properties, options);
    case 'Polygon':
        if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
        if (segmentIndex < 0) segmentIndex = coords[geometryIndex].length + segmentIndex - 1;
        return helpers.lineString([coords[geometryIndex][segmentIndex], coords[geometryIndex][segmentIndex + 1]], properties, options);
    case 'MultiLineString':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (segmentIndex < 0) segmentIndex = coords[multiFeatureIndex].length + segmentIndex - 1;
        return helpers.lineString([coords[multiFeatureIndex][segmentIndex], coords[multiFeatureIndex][segmentIndex + 1]], properties, options);
    case 'MultiPolygon':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (geometryIndex < 0) geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
        if (segmentIndex < 0) segmentIndex = coords[multiFeatureIndex][geometryIndex].length - segmentIndex - 1;
        return helpers.lineString([coords[multiFeatureIndex][geometryIndex][segmentIndex], coords[multiFeatureIndex][geometryIndex][segmentIndex + 1]], properties, options);
    }
    throw new Error('geojson is invalid');
}

/**
 * Finds a particular Point from a GeoJSON using `@turf/meta` indexes.
 *
 * Negative indexes are permitted.
 *
 * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {number} [options.featureIndex=0] Feature Index
 * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
 * @param {number} [options.geometryIndex=0] Geometry Index
 * @param {number} [options.coordIndex=0] Coord Index
 * @param {Object} [options.properties={}] Translate Properties to output Point
 * @param {BBox} [options.bbox={}] Translate BBox to output Point
 * @param {number|string} [options.id={}] Translate Id to output Point
 * @returns {Feature<Point>} 2-vertex GeoJSON Feature Point
 * @example
 * var multiLine = turf.multiLineString([
 *     [[10, 10], [50, 30], [30, 40]],
 *     [[-10, -10], [-50, -30], [-30, -40]]
 * ]);
 *
 * // First Segment (defaults are 0)
 * turf.findPoint(multiLine);
 * // => Feature<Point<[10, 10]>>
 *
 * // First Segment of the 2nd Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: 1});
 * // => Feature<Point<[-10, -10]>>
 *
 * // Last Segment of last Multi-Feature
 * turf.findPoint(multiLine, {multiFeatureIndex: -1, coordIndex: -1});
 * // => Feature<Point<[-30, -40]>>
 */
function findPoint(geojson, options) {
    // Optional Parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var featureIndex = options.featureIndex || 0;
    var multiFeatureIndex = options.multiFeatureIndex || 0;
    var geometryIndex = options.geometryIndex || 0;
    var coordIndex = options.coordIndex || 0;

    // Find FeatureIndex
    var properties = options.properties;
    var geometry;

    switch (geojson.type) {
    case 'FeatureCollection':
        if (featureIndex < 0) featureIndex = geojson.features.length + featureIndex;
        properties = properties || geojson.features[featureIndex].properties;
        geometry = geojson.features[featureIndex].geometry;
        break;
    case 'Feature':
        properties = properties || geojson.properties;
        geometry = geojson.geometry;
        break;
    case 'Point':
    case 'MultiPoint':
        return null;
    case 'LineString':
    case 'Polygon':
    case 'MultiLineString':
    case 'MultiPolygon':
        geometry = geojson;
        break;
    default:
        throw new Error('geojson is invalid');
    }

    // Find Coord Index
    if (geometry === null) return null;
    var coords = geometry.coordinates;
    switch (geometry.type) {
    case 'Point':
        return helpers.point(coords, properties, options);
    case 'MultiPoint':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        return helpers.point(coords[multiFeatureIndex], properties, options);
    case 'LineString':
        if (coordIndex < 0) coordIndex = coords.length + coordIndex;
        return helpers.point(coords[coordIndex], properties, options);
    case 'Polygon':
        if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
        if (coordIndex < 0) coordIndex = coords[geometryIndex].length + coordIndex;
        return helpers.point(coords[geometryIndex][coordIndex], properties, options);
    case 'MultiLineString':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (coordIndex < 0) coordIndex = coords[multiFeatureIndex].length + coordIndex;
        return helpers.point(coords[multiFeatureIndex][coordIndex], properties, options);
    case 'MultiPolygon':
        if (multiFeatureIndex < 0) multiFeatureIndex = coords.length + multiFeatureIndex;
        if (geometryIndex < 0) geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
        if (coordIndex < 0) coordIndex = coords[multiFeatureIndex][geometryIndex].length - coordIndex;
        return helpers.point(coords[multiFeatureIndex][geometryIndex][coordIndex], properties, options);
    }
    throw new Error('geojson is invalid');
}

exports.coordEach = coordEach;
exports.coordReduce = coordReduce;
exports.propEach = propEach;
exports.propReduce = propReduce;
exports.featureEach = featureEach;
exports.featureReduce = featureReduce;
exports.coordAll = coordAll;
exports.geomEach = geomEach;
exports.geomReduce = geomReduce;
exports.flattenEach = flattenEach;
exports.flattenReduce = flattenReduce;
exports.segmentEach = segmentEach;
exports.segmentReduce = segmentReduce;
exports.lineEach = lineEach;
exports.lineReduce = lineReduce;
exports.findSegment = findSegment;
exports.findPoint = findPoint;

},{"@turf/helpers":35}],35:[function(require,module,exports){
arguments[4][7][0].apply(exports,arguments)
},{"dup":7}],36:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
var bearing_1 = require("@turf/bearing");
var distance_1 = require("@turf/distance");
var destination_1 = require("@turf/destination");
var line_intersect_1 = require("@turf/line-intersect");
var meta_1 = require("@turf/meta");
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
/**
 * Takes a {@link Point} and a {@link LineString} and calculates the closest Point on the (Multi)LineString.
 *
 * @name nearestPointOnLine
 * @param {Geometry|Feature<LineString|MultiLineString>} lines lines to snap to
 * @param {Geometry|Feature<Point>|number[]} pt point to snap from
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
 * @returns {Feature<Point>} closest point on the `line` to `point`. The properties object will contain three values: `index`: closest point was found on nth line part, `dist`: distance between pt and the closest point, `location`: distance along the line between start and the closest point.
 * @example
 * var line = turf.lineString([
 *     [-77.031669, 38.878605],
 *     [-77.029609, 38.881946],
 *     [-77.020339, 38.884084],
 *     [-77.025661, 38.885821],
 *     [-77.021884, 38.889563],
 *     [-77.019824, 38.892368]
 * ]);
 * var pt = turf.point([-77.037076, 38.884017]);
 *
 * var snapped = turf.nearestPointOnLine(line, pt, {units: 'miles'});
 *
 * //addToMap
 * var addToMap = [line, pt, snapped];
 * snapped.properties['marker-color'] = '#00f';
 */
function nearestPointOnLine(lines, pt, options) {
    if (options === void 0) { options = {}; }
    var closestPt = helpers_1.point([Infinity, Infinity], {
        dist: Infinity
    });
    var length = 0.0;
    meta_1.flattenEach(lines, function (line) {
        var coords = invariant_1.getCoords(line);
        for (var i = 0; i < coords.length - 1; i++) {
            //start
            var start = helpers_1.point(coords[i]);
            start.properties.dist = distance_1.default(pt, start, options);
            //stop
            var stop_1 = helpers_1.point(coords[i + 1]);
            stop_1.properties.dist = distance_1.default(pt, stop_1, options);
            // sectionLength
            var sectionLength = distance_1.default(start, stop_1, options);
            //perpendicular
            var heightDistance = Math.max(start.properties.dist, stop_1.properties.dist);
            var direction = bearing_1.default(start, stop_1);
            var perpendicularPt1 = destination_1.default(pt, heightDistance, direction + 90, options);
            var perpendicularPt2 = destination_1.default(pt, heightDistance, direction - 90, options);
            var intersect = line_intersect_1.default(helpers_1.lineString([perpendicularPt1.geometry.coordinates, perpendicularPt2.geometry.coordinates]), helpers_1.lineString([start.geometry.coordinates, stop_1.geometry.coordinates]));
            var intersectPt = null;
            if (intersect.features.length > 0) {
                intersectPt = intersect.features[0];
                intersectPt.properties.dist = distance_1.default(pt, intersectPt, options);
                intersectPt.properties.location = length + distance_1.default(start, intersectPt, options);
            }
            if (start.properties.dist < closestPt.properties.dist) {
                closestPt = start;
                closestPt.properties.index = i;
                closestPt.properties.location = length;
            }
            if (stop_1.properties.dist < closestPt.properties.dist) {
                closestPt = stop_1;
                closestPt.properties.index = i + 1;
                closestPt.properties.location = length + sectionLength;
            }
            if (intersectPt && intersectPt.properties.dist < closestPt.properties.dist) {
                closestPt = intersectPt;
                closestPt.properties.index = i;
            }
            // update length
            length += sectionLength;
        }
    });
    return closestPt;
}
exports.default = nearestPointOnLine;

},{"@turf/bearing":4,"@turf/destination":13,"@turf/distance":14,"@turf/helpers":15,"@turf/invariant":17,"@turf/line-intersect":20,"@turf/meta":37}],37:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],38:[function(require,module,exports){
"use strict";
var __importDefault = (this && this.__importDefault) || function (mod) {
    return (mod && mod.__esModule) ? mod : { "default": mod };
}
Object.defineProperty(exports, "__esModule", { value: true });
// Taken from http://geomalgorithms.com/a02-_lines.html
var distance_1 = __importDefault(require("@turf/distance"));
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
var meta_1 = require("@turf/meta");
var rhumb_distance_1 = __importDefault(require("@turf/rhumb-distance"));
/**
 * Returns the minimum distance between a {@link Point} and a {@link LineString}, being the distance from a line the
 * minimum distance between the point and any segment of the `LineString`.
 *
 * @name pointToLineDistance
 * @param {Feature<Point>|Array<number>} pt Feature or Geometry
 * @param {Feature<LineString>} line GeoJSON Feature or Geometry
 * @param {Object} [options={}] Optional parameters
 * @param {string} [options.units="kilometers"] can be anything supported by turf/convertLength
 * (ex: degrees, radians, miles, or kilometers)
 * @param {string} [options.method="geodesic"] wether to calculate the distance based on geodesic (spheroid) or
 * planar (flat) method. Valid options are 'geodesic' or 'planar'.
 * @returns {number} distance between point and line
 * @example
 * var pt = turf.point([0, 0]);
 * var line = turf.lineString([[1, 1],[-1, 1]]);
 *
 * var distance = turf.pointToLineDistance(pt, line, {units: 'miles'});
 * //=69.11854715938406
 */
function pointToLineDistance(pt, line, options) {
    if (options === void 0) { options = {}; }
    // Optional parameters
    if (!options.method) {
        options.method = "geodesic";
    }
    if (!options.units) {
        options.units = "kilometers";
    }
    // validation
    if (!pt) {
        throw new Error("pt is required");
    }
    if (Array.isArray(pt)) {
        pt = helpers_1.point(pt);
    }
    else if (pt.type === "Point") {
        pt = helpers_1.feature(pt);
    }
    else {
        invariant_1.featureOf(pt, "Point", "point");
    }
    if (!line) {
        throw new Error("line is required");
    }
    if (Array.isArray(line)) {
        line = helpers_1.lineString(line);
    }
    else if (line.type === "LineString") {
        line = helpers_1.feature(line);
    }
    else {
        invariant_1.featureOf(line, "LineString", "line");
    }
    var distance = Infinity;
    var p = pt.geometry.coordinates;
    meta_1.segmentEach(line, function (segment) {
        var a = segment.geometry.coordinates[0];
        var b = segment.geometry.coordinates[1];
        var d = distanceToSegment(p, a, b, options);
        if (d < distance) {
            distance = d;
        }
    });
    return helpers_1.convertLength(distance, "degrees", options.units);
}
/**
 * Returns the distance between a point P on a segment AB.
 *
 * @private
 * @param {Array<number>} p external point
 * @param {Array<number>} a first segment point
 * @param {Array<number>} b second segment point
 * @param {Object} [options={}] Optional parameters
 * @returns {number} distance
 */
function distanceToSegment(p, a, b, options) {
    var v = [b[0] - a[0], b[1] - a[1]];
    var w = [p[0] - a[0], p[1] - a[1]];
    var c1 = dot(w, v);
    if (c1 <= 0) {
        return calcDistance(p, a, { method: options.method, units: "degrees" });
    }
    var c2 = dot(v, v);
    if (c2 <= c1) {
        return calcDistance(p, b, { method: options.method, units: "degrees" });
    }
    var b2 = c1 / c2;
    var Pb = [a[0] + (b2 * v[0]), a[1] + (b2 * v[1])];
    return calcDistance(p, Pb, { method: options.method, units: "degrees" });
}
function dot(u, v) {
    return (u[0] * v[0] + u[1] * v[1]);
}
function calcDistance(a, b, options) {
    return options.method === "planar" ? rhumb_distance_1.default(a, b, options) : distance_1.default(a, b, options);
}
exports.default = pointToLineDistance;

},{"@turf/distance":14,"@turf/helpers":15,"@turf/invariant":17,"@turf/meta":39,"@turf/rhumb-distance":42}],39:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],40:[function(require,module,exports){
'use strict';

Object.defineProperty(exports, '__esModule', { value: true });

function _interopDefault (ex) { return (ex && (typeof ex === 'object') && 'default' in ex) ? ex['default'] : ex; }

var meta = require('@turf/meta');
var helpers = require('@turf/helpers');
var clone = _interopDefault(require('@turf/clone'));

/**
 * Converts a WGS84 GeoJSON object into Mercator (EPSG:900913) projection
 *
 * @name toMercator
 * @param {GeoJSON|Position} geojson WGS84 GeoJSON object
 * @param {Object} [options] Optional parameters
 * @param {boolean} [options.mutate=false] allows GeoJSON input to be mutated (significant performance increase if true)
 * @returns {GeoJSON} true/false
 * @example
 * var pt = turf.point([-71,41]);
 * var converted = turf.toMercator(pt);
 *
 * //addToMap
 * var addToMap = [pt, converted];
 */
function toMercator(geojson, options) {
    return convert(geojson, 'mercator', options);
}

/**
 * Converts a Mercator (EPSG:900913) GeoJSON object into WGS84 projection
 *
 * @name toWgs84
 * @param {GeoJSON|Position} geojson Mercator GeoJSON object
 * @param {Object} [options] Optional parameters
 * @param {boolean} [options.mutate=false] allows GeoJSON input to be mutated (significant performance increase if true)
 * @returns {GeoJSON} true/false
 * @example
 * var pt = turf.point([-7903683.846322424, 5012341.663847514]);
 * var converted = turf.toWgs84(pt);
 *
 * //addToMap
 * var addToMap = [pt, converted];
 */
function toWgs84(geojson, options) {
    return convert(geojson, 'wgs84', options);
}


/**
 * Converts a GeoJSON coordinates to the defined `projection`
 *
 * @private
 * @param {GeoJSON} geojson GeoJSON Feature or Geometry
 * @param {string} projection defines the projection system to convert the coordinates to
 * @param {Object} [options] Optional parameters
 * @param {boolean} [options.mutate=false] allows GeoJSON input to be mutated (significant performance increase if true)
 * @returns {GeoJSON} true/false
 */
function convert(geojson, projection, options) {
    // Optional parameters
    options = options || {};
    if (!helpers.isObject(options)) throw new Error('options is invalid');
    var mutate = options.mutate;

    // Validation
    if (!geojson) throw new Error('geojson is required');

    // Handle Position
    if (Array.isArray(geojson) && helpers.isNumber(geojson[0])) geojson = (projection === 'mercator') ? convertToMercator(geojson) : convertToWgs84(geojson);

    // Handle GeoJSON
    else {
        // Handle possible data mutation
        if (mutate !== true) geojson = clone(geojson);

        meta.coordEach(geojson, function (coord) {
            var newCoord = (projection === 'mercator') ? convertToMercator(coord) : convertToWgs84(coord);
            coord[0] = newCoord[0];
            coord[1] = newCoord[1];
        });
    }
    return geojson;
}

/**
 * Convert lon/lat values to 900913 x/y.
 * (from https://github.com/mapbox/sphericalmercator)
 *
 * @private
 * @param {Array<number>} lonLat WGS84 point
 * @returns {Array<number>} Mercator [x, y] point
 */
function convertToMercator(lonLat) {
    var D2R = Math.PI / 180,
        // 900913 properties
        A = 6378137.0,
        MAXEXTENT = 20037508.342789244;

    // compensate longitudes passing the 180th meridian
    // from https://github.com/proj4js/proj4js/blob/master/lib/common/adjust_lon.js
    var adjusted = (Math.abs(lonLat[0]) <= 180) ? lonLat[0] : (lonLat[0] - (sign(lonLat[0]) * 360));
    var xy = [
        A * adjusted * D2R,
        A * Math.log(Math.tan((Math.PI * 0.25) + (0.5 * lonLat[1] * D2R)))
    ];

    // if xy value is beyond maxextent (e.g. poles), return maxextent
    if (xy[0] > MAXEXTENT) xy[0] = MAXEXTENT;
    if (xy[0] < -MAXEXTENT) xy[0] = -MAXEXTENT;
    if (xy[1] > MAXEXTENT) xy[1] = MAXEXTENT;
    if (xy[1] < -MAXEXTENT) xy[1] = -MAXEXTENT;

    return xy;
}

/**
 * Convert 900913 x/y values to lon/lat.
 * (from https://github.com/mapbox/sphericalmercator)
 *
 * @private
 * @param {Array<number>} xy Mercator [x, y] point
 * @returns {Array<number>} WGS84 [lon, lat] point
 */
function convertToWgs84(xy) {
    // 900913 properties.
    var R2D = 180 / Math.PI;
    var A = 6378137.0;

    return [
        (xy[0] * R2D / A),
        ((Math.PI * 0.5) - 2.0 * Math.atan(Math.exp(-xy[1] / A))) * R2D
    ];
}

/**
 * Returns the sign of the input, or zero
 *
 * @private
 * @param {number} x input
 * @returns {number} -1|0|1 output
 */
function sign(x) {
    return (x < 0) ? -1 : (x > 0) ? 1 : 0;
}

exports.toMercator = toMercator;
exports.toWgs84 = toWgs84;

},{"@turf/clone":12,"@turf/helpers":41,"@turf/meta":34}],41:[function(require,module,exports){
arguments[4][7][0].apply(exports,arguments)
},{"dup":7}],42:[function(require,module,exports){
"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
// https://en.wikipedia.org/wiki/Rhumb_line
var helpers_1 = require("@turf/helpers");
var invariant_1 = require("@turf/invariant");
/**
 * Calculates the distance along a rhumb line between two {@link Point|points} in degrees, radians,
 * miles, or kilometers.
 *
 * @name rhumbDistance
 * @param {Coord} from origin point
 * @param {Coord} to destination point
 * @param {Object} [options] Optional parameters
 * @param {string} [options.units="kilometers"] can be degrees, radians, miles, or kilometers
 * @returns {number} distance between the two points
 * @example
 * var from = turf.point([-75.343, 39.984]);
 * var to = turf.point([-75.534, 39.123]);
 * var options = {units: 'miles'};
 *
 * var distance = turf.rhumbDistance(from, to, options);
 *
 * //addToMap
 * var addToMap = [from, to];
 * from.properties.distance = distance;
 * to.properties.distance = distance;
 */
function rhumbDistance(from, to, options) {
    if (options === void 0) { options = {}; }
    var origin = invariant_1.getCoord(from);
    var destination = invariant_1.getCoord(to);
    // compensate the crossing of the 180th meridian (https://macwright.org/2016/09/26/the-180th-meridian.html)
    // solution from https://github.com/mapbox/mapbox-gl-js/issues/3250#issuecomment-294887678
    destination[0] += (destination[0] - origin[0] > 180) ? -360 : (origin[0] - destination[0] > 180) ? 360 : 0;
    var distanceInMeters = calculateRhumbDistance(origin, destination);
    var distance = helpers_1.convertLength(distanceInMeters, "meters", options.units);
    return distance;
}
/**
 * Returns the distance travelling from ‘this’ point to destination point along a rhumb line.
 * Adapted from Geodesy: https://github.com/chrisveness/geodesy/blob/master/latlon-spherical.js
 *
 * @private
 * @param   {Array<number>} origin point.
 * @param   {Array<number>} destination point.
 * @param   {number} [radius=6371e3] - (Mean) radius of earth (defaults to radius in metres).
 * @returns {number} Distance in km between this point and destination point (same units as radius).
 *
 * @example
 *     var p1 = new LatLon(51.127, 1.338);
 *     var p2 = new LatLon(50.964, 1.853);
 *     var d = p1.distanceTo(p2); // 40.31 km
 */
function calculateRhumbDistance(origin, destination, radius) {
    // φ => phi
    // λ => lambda
    // ψ => psi
    // Δ => Delta
    // δ => delta
    // θ => theta
    radius = (radius === undefined) ? helpers_1.earthRadius : Number(radius);
    // see www.edwilliams.org/avform.htm#Rhumb
    var R = radius;
    var phi1 = origin[1] * Math.PI / 180;
    var phi2 = destination[1] * Math.PI / 180;
    var DeltaPhi = phi2 - phi1;
    var DeltaLambda = Math.abs(destination[0] - origin[0]) * Math.PI / 180;
    // if dLon over 180° take shorter rhumb line across the anti-meridian:
    if (DeltaLambda > Math.PI) {
        DeltaLambda -= 2 * Math.PI;
    }
    // on Mercator projection, longitude distances shrink by latitude; q is the 'stretch factor'
    // q becomes ill-conditioned along E-W line (0/0); use empirical tolerance to avoid it
    var DeltaPsi = Math.log(Math.tan(phi2 / 2 + Math.PI / 4) / Math.tan(phi1 / 2 + Math.PI / 4));
    var q = Math.abs(DeltaPsi) > 10e-12 ? DeltaPhi / DeltaPsi : Math.cos(phi1);
    // distance is pythagoras on 'stretched' Mercator projection
    var delta = Math.sqrt(DeltaPhi * DeltaPhi + q * q * DeltaLambda * DeltaLambda); // angular distance in radians
    var dist = delta * R;
    return dist;
}
exports.default = rhumbDistance;

},{"@turf/helpers":15,"@turf/invariant":17}],43:[function(require,module,exports){
/* The Agentmap class, which turns a Leaflet map into a simulation platform. */

let lineSlice = require('@turf/line-slice').default,
length = require('@turf/length').default;

/**
 * The main class for building, storing, simulating, and manipulating agent-based models on Leaflet maps.
 *
 * @class Agentmap
 * @param {object} map - A Leaflet Map instance.
 * @param {number} [animation_interval=1] - The number of steps agents must move before being redrawn. Given 1, they will be redrawn after every step. Given 0, the animation will not update at all. 1 by default. Must be a nonnegative integer.
 * @property {object} map - A Leaflet Map instance.
 * @property {FeatureGroup} agents - A featureGroup containing all agents.
 * @property {FeatureGroup} units - A featureGroup containing all units.
 * @property {FeatureGroup} streets - A featureGroup containing all streets.
 * @property {object} state - Properties detailing the state of the simulation process.
 * @property {boolean} state.running - Whether the simulation is running or not.
 * @property {boolean} state.paused - Whether the simulation is paused.
 * @property {?number} state.animation_frame_id - The id of the agentmap's update function in the queue of functions to call for the coming animation frame.
 * @property {?number} state.ticks - The number of ticks elapsed since the start of the simulation.
 * @property {number} animation_interval - The number of steps agents must move before being redrawn. Given 1, they will be redrawn after every step. Given 0, the animation will not update at all. 1 by default. Will be a nonnegative integer.
 * @property {?function} controller - User-defined function to be called on each update.
 */
Agentmap = function (map, animation_interval = 1) {
	Agentmap.checkAnimIntervalOption(animation_interval);

	this.map = map,
	this.units = null,
	this.streets = null,
	this.agents = null, 
	this.pathfinder = null,
	this.state = {
		running: false,
		paused: false,
		animation_frame_id: null,
		ticks: null,
	},
	this.controller = function() {},
	this.animation_interval = animation_interval
};

/**
 * Change the animation interval of the simulation & redraw the agents.
 *
 * @param {number} animation_interval - The desired animation interval to give the simulation. Must be a nonnegative integer.
 */
Agentmap.prototype.setAnimationInterval = function(animation_interval) {
	Agentmap.checkAnimIntervalOption(animation_interval);

	this.animation_interval = animation_interval;

	this.agents.eachLayer(agent => agent.setLatLng(agent._latlng));
}

/**
 * Check whether the animation interval option provided is valid.
 * @private
 *
 * @param {number} animation_interval - An input specifying an animation interval distance.
 */
Agentmap.checkAnimIntervalOption = function(animation_interval) {
	if (!Number.isInteger(animation_interval) && animation_interval >= 0) {
		throw new Error("The animation_interval must be a non-negative integer!");
	}
}

/**
 * Get an animation frame, have the agents update & get ready to be drawn, and keep doing that until paused or reset.
 */
Agentmap.prototype.run = function() {
	if (this.state.running === false) {
		this.state.running = true;

		let animation_update = (function (rAF_time) {
			
			if (this.state.paused === true) {
				this.state.paused = false;
			}
			this.state.animation_frame_id = L.Util.requestAnimFrame(animation_update);
			this.update();
		}).bind(this);

		this.state.animation_frame_id = L.Util.requestAnimFrame(animation_update);
	}
}

/**
 * Update the simulation at the given time.
 * @private
 */
Agentmap.prototype.update = function() {
	if (this.state.ticks === null) {
		this.state.ticks = 0;
	}

	//Execute user-provided per-tick instructions for the agentmap.
	this.controller();

	//Execute user-provided per-tick instructions for each agent.
	this.agents.eachLayer(function(agent) {
		agent.controller();
	});
	
	this.state.ticks += 1;
};

/**
* Stop the animation, reset the animation state properties, and delete the features.
*/
Agentmap.prototype.clear = function() {
	L.Util.cancelAnimFrame(this.state.animation_frame_id);
	this.state.running = false,
	this.state.paused = false,
	this.state.animation_frame_id = null,
	this.state.ticks = null,
	
	this.agents.clearLayers();
	this.streets.clearLayers();
	this.units.clearLayers();
};

/** 
 * Stop the animation, stop updating the agents.
 */
Agentmap.prototype.pause = function() {
	L.Util.cancelAnimFrame(this.state.animation_frame_id);
	this.state.running = false,
	this.state.paused = true;
};

/**
 * Get a point through which an agent can exit/enter a unit.
 *
 * @param {number} unit_id - The unique ID of the unit whose door you want.
 * @returns {LatLng} - The coordinates of the center point of the segment of the unit parallel to the street.
 */
Agentmap.prototype.getUnitDoor = function(unit_id) {
	let unit = this.units.getLayer(unit_id);
	
	if (typeof(unit) === "undefined") {
		throw new Error("No unit with the specified ID exists.");
	}
	
	let unit_spec = unit.getLatLngs()[0],
	corner_a = unit_spec[0],
	corner_b = unit_spec[1],
	door = 	L.latLngBounds(corner_a, corner_b).getCenter();
	
	return door;
};

/**
 * Get the point on the adjacent street in front of the unit's door.
 *
 * @param {number} unit_id - The unique ID of the unit whose door's corresponding point on the street you want.
 * @returns {LatLng} - The coordinates point of the adjacent street directly in front of unit's door.
 */
Agentmap.prototype.getStreetNearDoor = function(unit_id) {
	let unit = this.units.getLayer(unit_id);
	
	if (typeof(unit) === "undefined") {
		throw new Error("No unit with the specified ID exists.");
	}
	
	let unit_anchors = L.A.reversedCoordinates(unit.street_anchors),
	street_point = L.latLngBounds(...unit_anchors).getCenter();
	
	return street_point;
};

/**
 * Given a unit and a pair of coordinates between 0 and 1, return a corresponding point inside the unit, offset from its first corner along the street.
 * 
 * @param {number} unit_id - The unique ID of the unit whose interior point you want.
 * @param {number} x - A point between 0 and 1 representing a position along the width of a unit.
 * @param {number} y - A point between 0 and 1 representing a position along the depth of a unit.
 * @returns {LatLng} - The global coordinates of the specified position within the unit.
 */
Agentmap.prototype.getUnitPoint = function(unit_id, x, y) {
	if (x < 0 || x > 1 || y < 0 || y > 1) {
		throw new Error("x and y must both be between 0 and 1!");
	}
	
	let unit = this.units.getLayer(unit_id),
	unit_corners = unit.getLatLngs()[0],
	front_right = unit_corners[0],
	front_left = unit_corners[1],
	back_right = unit_corners[3],
	front_length = front_left.lng - front_right.lng,
	side_length = back_right.lng - front_right.lng,
	front_slope = (front_right.lat - front_left.lat) / (front_right.lng - front_left.lng),
	side_slope = (front_right.lat - back_right.lat) / (front_right.lng - back_right.lng);
	
	//Get the coordinate of the position along the front (x) axis.
	let lng_along_front = front_right.lng + front_length * x,
	lat_along_front = front_right.lat + (front_length * x) * front_slope,
	point_along_front = L.latLng(lat_along_front, lng_along_front);
	
	//From the position on the front axis, get the coordinate of a position along a line perpendicular to the front and 
	//parallel to the side (y) axis.
	let lng_along_side = point_along_front.lng + side_length * y,
	lat_along_side = point_along_front.lat + (side_length * y) * side_slope,
	point_in_depth = L.latLng(lat_along_side, lng_along_side);

	return point_in_depth;
}

/**
 * Given a point on a street, find the nearest intersection on that street (with any other street).
 * 
 * @param {LatLng} lat_lng - The coordinates of the point on the street to search from.
 * @param {Place} place - A place object corresponding to the street.
 * @returns {LatLng} - The coordinates of the nearest intersection.
 */
Agentmap.prototype.getNearestIntersection = function(lat_lng, place) {
	let street_id,
	street_feature;

	if (place.type === "street") {
		street_id = place.id;
	}
	else {
		throw new Error("place must be a street!");
	}

	street_feature = this.streets.getLayer(street_id).toGeoJSON();
		
	let intersections = this.streets.getLayer(street_id).intersections,
	intersection_points = [],
	intersection_distances = [];

	for (let intersection in intersections) { 
		for (let cross_point of intersections[intersection]) {
			let intersection_point = cross_point[0],
			distance = lat_lng.distanceTo(intersection_point);

			/* More precise, but slower, distance detection -- not necessary yet. 
			 	let start_coords = L.A.pointToCoordinateArray(lat_lng);
				intersection_coords = L.A.pointToCoordinateArray(intersection_point),
				segment = lineSlice(start_coords, intersection_coords, street_feature),
				distance = length(segment); 
			*/
			
			intersection_points.push(intersection_point);
			intersection_distances.push(distance);
		}
	}
	
	let smallest_distance = Math.min(...intersection_distances),
	smallest_distance_index = intersection_distances.indexOf(smallest_distance),
	closest_intersection_point = L.latLng(intersection_points[smallest_distance_index]);
	
	return closest_intersection_point;
}

/**
 * Since units may take a noticeably long time to generate while typically staying the same over simulations,
 * downloadUnits makes it easy to get a JS file containing the units object, so it can be included with an
 * AgentMaps app and imported into Agentmap.buildingify so they will not need to be regenerated.
 */
Agentmap.prototype.downloadUnits = function() {
	let file_content = "let units_data = ",
	units_json = this.units.toGeoJSON(20);
	file_content += JSON.stringify(units_json),
	file = new Blob([file_content]);

	var element = document.createElement("a");
	element.setAttribute("href", URL.createObjectURL(file)),
	element.setAttribute("download", "units_data.js"),
	element.style.display = "none";
	document.body.appendChild(element);
	
	element.click();
	
	document.body.removeChild(element);
}

/**
 * Generates an agentmap for the given map.
 *
 * @name agentmap
 * @param {object} map - A Leaflet Map instance.
 * @returns {object} - An Agentmap instance.
 */
function agentmapFactory(map) {
	return new Agentmap(map);
}

/**
 * Returns the number of layers in a Leaflet layer group.
 *
 * @memberof L.LayerGroup
 */
function layerCount() {
	return this.getLayers().length;
}

L.LayerGroup.include({count: layerCount});

exports.Agentmap = Agentmap,
exports.agentmap = agentmapFactory;

},{"@turf/length":18,"@turf/line-slice":24}],44:[function(require,module,exports){
/* Here we define agentify, the agent base class, and everything they uniquely rely on. */

let centroid = require('@turf/centroid').default,
buffer = require('@turf/buffer').default,
booleanPointInPolygon = require('@turf/boolean-point-in-polygon').default,
along = require('@turf/along').default,
nearestPointOnLine = require('@turf/nearest-point-on-line').default,
lineSlice = require('@turf/line-slice').default,
length = require('@turf/length').default,
lineString = require('@turf/helpers').lineString,
bearing = require('@turf/bearing').default,
destination = require('@turf/destination').default,
Agentmap = require('./agentmap').Agentmap,
encodeLatLng = require('./routing').encodeLatLng;

/**
 * The main class representing individual agents, using Leaflet class system.
 * @private
 *
 * @class Agent
 */
let Agent = {};

/**
 * Constructor for the Agent class, using Leaflet class system.
 *
 * @name Agent
 * @constructor
 * @param {LatLng} lat_lng - A pair of coordinates to place the agent at.
 * @param {Object} options - An array of options for the agent, namely its layer.
 * @param {Agentmap} agentmap - The agentmap instance in which the agent exists.
 * @property {number} feature.AgentMap_id - The agent's instance id, so it can be accessed from inside the Leaflet layer. To avoid putting the actual instance inside the feature object.
 * @property {Agentmap} agentmap - The agentmap instance in which the agent exists.
 * @property {Place} place - A place object specifying where the agent is currently at.
 * @property {number} [steps_made=0] - The number of steps the agent has moved since the beginning.
 * @property {Object} this.trip - Properties detailing information about the agent's trip that change sometimes, but needs to be accessed by future updates.
 * @property {boolean} this.trip.moving - Whether the agent currently moving.
 * @property {boolean} this.trip.paused - Whether the agent should be allowed to move along its trip.
 * @property {?Point} this.trip.current_point - The point where the agent is currently located.
 * @property {?Point} this.trip.goal_point - The point where the agent is traveling to.
 * @property {?number} this.trip.lat_dir - The latitudinal direction. -1 if traveling to lower latitude (down), 1 if traveling to higher latitude (up).
 * @property {?number} this.trip.lng_dir - The longitudinal direction. -1 if traveling to lesser longitude (left), 1 if traveling to greater longitude (right).
 * @property {?number} this.trip.speed - The speed that the agent should travel, in meters per tick.
 * @property {?number} this.trip.angle - The angle between the current point and the goal.
 * @property {?number} this.trip.slope - The slope of the line segment formed by the two points between which the agent is traveling at this time during its trip.
 * @property {Array} this.trip.path - A sequence of LatLngs; the agent will move from one to the next, popping each one off after it arrives until the end of the street; or, until the trip is changed/reset.
 * @property {?function} controller - User-defined function to be called on each update (each tick).
 * @property {?function} fine_controller - User-defined function to be called before & after each movemnt (on each step an agent performs during a tick).
 */
Agent.initialize = function(lat_lng, options, agentmap) {
	this.agentmap = agentmap,
	this.place = null,
	this.steps_made = 0,
	this.trip = {
		paused: false,
		moving: false,
		current_point: null,
		goal_point: null,
		lat_dir: null,
		lng_dir: null,
		slope: null,
		angle: null,
		speed: null,
		path: [],
	},
	this.controller = function() {},
	this.fine_controller = function() {};

	L.CircleMarker.prototype.initialize.call(this, lat_lng, options);
}

/**
 * Reset all the properties of its trip, but don't change whether it's allowed to be traveling or not.
 * @memberof Agent
 * @instance
 */
Agent.resetTrip = function() {
	for (let key in this.trip) {
		this.trip[key] =
			key === "paused" ? false :
			key === "moving" ? false :
			key === "path" ? [] :
			null;
	}
};

/**
 * Set the agent up to start traveling along the path specified in the agent's trip..
 * @memberof Agent
 * @instance
 */
Agent.startTrip = function() {
	if (this.trip.path.length > 0) {
		this.travelTo(this.trip.path[0]);
	}
};

/**
 * Stop the agent where it is along its trip.
 * @memberof Agent
 * @instance
 */
Agent.pauseTrip = function() {
	this.trip.paused = true;
};

/**
 * Have the agent continue from where it was left off along its trip.
 * @memberof Agent
 * @instance
 */
Agent.resumeTrip = function() {
	this.trip.paused = false;
};

/**
 * Set the agent to travel to some point on the map.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {LatLng} goal_point - The point to which the agent should travel.
 */
Agent.travelTo = function(goal_point) {
	this.trip.current_point = this.getLatLng(),
	this.trip.goal_point = goal_point,

	//Negating so that neg result corresponds to the goal being rightward/above, pos result to it being leftward/below.
	this.trip.lat_dir = Math.sign(- (this.trip.current_point.lat - this.trip.goal_point.lat)),
	this.trip.lng_dir = Math.sign(- (this.trip.current_point.lng - this.trip.goal_point.lng)),

	this.trip.angle = bearing(L.A.pointToCoordinateArray(this.trip.current_point), L.A.pointToCoordinateArray(this.trip.goal_point));
	this.trip.slope = Math.abs((this.trip.current_point.lat - this.trip.goal_point.lat) / (this.trip.current_point.lng - this.trip.goal_point.lng));
	this.trip.speed = this.trip.goal_point.speed;

	//If the agent won't be at any particular place at least until it reaches its next goal, mark its place as unanchored.
	if (this.trip.path[0].new_place.type === "unanchored" || this.trip.path[0].move_directly === true) {
		this.place = {type: "unanchored"};
	}
};

/**
 * Given the agent's currently scheduledthis.trips (its path), get the place from which a newthis.trip should start (namely, the end of the current path).
 * That is: If there's already a path in queue, start the new path from the end of the existing one.
 * @memberof Agent
 * @instance
 * @private
 *
 * @returns {Place} - The place where a newthis.trip should start.
 */
Agent.newTripStartPlace = function() {
	if (this.trip.path.length === 0) {
		start_place = this.place;
	}
	else {
		start_place = this.trip.path[this.trip.path.length - 1].new_place;
	}

	return start_place;
}

/**
 * Schedule the agent to travel to a point within the unit he is in.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {LatLng} goal_lat_lng - LatLng coordinate object for a point in the same unit the agent is in.
 * @param {number} speed - The speed that the agent should travel, in meters per tick.
 */
Agent.setTravelInUnit = function(goal_lat_lng, goal_place, speed) {
	goal_lat_lng.new_place = goal_place,
	goal_lat_lng.speed = speed;
	this.trip.path.push(goal_lat_lng);
};

/**
 * Schedule the agent to travel directly from any point (e.g. of a street or unit) to a point (e.g. of another street or unit).
 * @memberof Agent
 * @instance
 *
 * @param {LatLng} goal_lat_lng - The point within the place to which the agent is to travel.
 * @param {Place} goal_place - The place to which the agent will travel.
 * @param {number} [speed=1] - The speed in meters per tick that the agent should try to travel. Must be >= .1.
 * @param {Boolean} [move_directly=false] - Whether to ignore the streets & roads and move directly to the goal.
 * @param {Boolean} [replace_trip=false] - Whether to empty the currently scheduled path and replace it with this new trip; false by default (the new trip is
 * simply appended to the current scheduled path).
 */
Agent.setTravelToPlace = function(goal_lat_lng, goal_place, speed = 1, move_directly = false, replace_trip = false) {
	this.checkSpeed(speed);

	let start_place = this.newTripStartPlace();
	goal_lat_lng = L.latLng(goal_lat_lng);

	if (replace_trip === true) {
		this.resetTrip();
	}

	//If either the agent is already unanchored or its goal is unanchored, just schedule it to move directly to its goal.
	if (start_place.type === "unanchored" || goal_place.type === "unanchored" || move_directly === true) {
		// if(!move_directly){debugger;}
		let goal = goal_lat_lng;
		goal.new_place = goal_place,
		goal.move_directly = true,
		goal.speed = speed;

		this.trip.path.push(goal);

		return;
	}

	let goal_layer = this.agentmap.units.getLayer(goal_place.id) || this.agentmap.streets.getLayer(goal_place.id);

	//If the goal isn't unanchored, see if it's a street or a unit and schedule the agent appropriately.
	if (goal_layer) {
		let goal_coords = L.A.pointToCoordinateArray(goal_lat_lng);

		//Buffering so that points on the perimeter, like the door, are captured.
		//Also expands street lines into thin polygons (booleanPointInPolygon requires polys).
		//Might be more efficient to generate the door so that it's slightly inside the area.
		let goal_polygon = buffer(goal_layer.toGeoJSON(), .001);

		if (booleanPointInPolygon(goal_coords, goal_polygon)) {
			if (start_place.type === "unit" && goal_place.type === "unit" && start_place.id === goal_place.id) {
				this.setTravelInUnit(goal_lat_lng, goal_place, speed);
				return;
			}
			//Move to the street if it's starting at a unit and its goal is elsewhere.
			else if (start_place.type === "unit") {
				let start_unit_door = this.agentmap.getUnitDoor(start_place.id);
				start_unit_door.new_place = start_place,
				start_unit_door.speed = speed;
				this.trip.path.push(start_unit_door);

				let start_unit_street_id = this.agentmap.units.getLayer(start_place.id).street_id,
				start_unit_street_point = this.agentmap.getStreetNearDoor(start_place.id);
				start_unit_street_point.new_place = { type: "street", id: start_unit_street_id },
				start_unit_street_point.speed = speed;
				this.trip.path.push(start_unit_street_point);
			}

			if (goal_place.type === "unit") {
				let goal_street_point = this.agentmap.getStreetNearDoor(goal_place.id),
				goal_street_point_place = { type: "street", id: this.agentmap.units.getLayer(goal_place.id).street_id };

				//Move to the point on the street closest to the goal unit...
				this.setTravelAlongStreet(goal_street_point, goal_street_point_place, speed);

				//Move from that point into the unit.
				let goal_door = this.agentmap.getUnitDoor(goal_place.id);
				goal_door.new_place = goal_place,
				goal_door.speed = speed;
				this.trip.path.push(goal_door)
				this.setTravelInUnit(goal_lat_lng, goal_place, speed);
			}
			else if (goal_place.street === "number") {
				this.setTravelAlongStreet(goal_lat_lng, goal_place, speed);
			}
		}
		else {
			throw new Error("The goal_lat_lng is not inside of the polygon of the goal_place!");
		}
	}
	else {
		throw new Error("No place exists matching the specified goal_place!");
	}
};

Agent.scheduleTrip = Agent.setTravelToPlace;

/**
 * Schedule the agent to travel to a point along the streets, via streets.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {LatLng} goal_lat_lng - The coordinates of a point on a street to which the agent should travel.
 * @param {Place} goal_place - The place to which the agent will travel. Must be a street.
 * @param {number} speed - The speed that the agent should travel, in meters per tick.
 */
Agent.setTravelAlongStreet = function(goal_lat_lng, goal_place, speed) {
	let goal_coords,
	goal_street_id,
	goal_street_point,
	goal_street_feature,
	start_place = this.newTripStartPlace(),
	start_street_id,
	start_street_point,
	start_street_feature;

	if (start_place.type === "street" && goal_place.type === "street") {
		start_street_id = start_place.id,
		start_street_point = this.trip.path.length !== 0 ?
			this.trip.path[this.trip.path.length - 1] :
			this.getLatLng();
		start_street_point.new_place = {type: "street", id: start_street_id};

		goal_street_id = goal_place.id,
		goal_street_feature = this.agentmap.streets.getLayer(goal_street_id).feature,
		goal_coords = L.A.pointToCoordinateArray(goal_lat_lng),
		goal_street_point = L.latLng(nearestPointOnLine(goal_street_feature, goal_coords).geometry.coordinates.reverse());
		goal_street_point.new_place = goal_place;
	}
	else {
		throw new Error("Both the start and end places must be streets!");
	}

	if (start_street_id === goal_street_id) {
		this.setTravelOnSameStreet(start_street_point, goal_street_point, goal_street_feature, goal_street_id, speed);
	}
	//If the start and end points are on different streets, move from the start to its nearest intersection, then from there
	//to the intersection nearest to the end, and finally to the end.
	else {
		let start_nearest_intersection = this.agentmap.getNearestIntersection(start_street_point, start_place),
		goal_nearest_intersection = this.agentmap.getNearestIntersection(goal_street_point, goal_place);

		start_street_feature = this.agentmap.streets.getLayer(start_street_id).feature;

		this.setTravelOnStreetNetwork(start_street_point, goal_street_point, start_nearest_intersection, goal_nearest_intersection, speed);
	}
};

/**
 * Schedule the agent to travel between two points on the same street.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param start_lat_lng {LatLng} - The coordinates of the point on the street from which the agent will be traveling.
 * @param goal_lat_lng {LatLng} - The coordinates of the point on the street to which the agent should travel.
 * @param street_feature {Feature} - A GeoJSON object representing an OpenStreetMap street.
 * @param street_id {number} - The ID of the street in the streets layerGroup.
 * @param {number} speed - The speed that the agent should travel, in meters per tick.
 */
Agent.setTravelOnSameStreet = function(start_lat_lng, goal_lat_lng, street_feature, street_id, speed) {
	//lineSlice, regardless of the specified starting point, will give a segment with the same coordinate order
	//as the original lineString array. So, if the goal point comes earlier in the array (e.g. it's on the far left),
	//it'll end up being the first point in the path, instead of the last, and the agent will move to it directly,
	//ignoring the street points that should come before it. It would then travel along the street from the goal point
	//to its original point (backwards).
	//To fix this, I'm reversing the order of the coordinates in the segment if the last point in the line is closer
	//to the agent's starting point than the first point on the line (implying the last point in the array is the starting
	//point, not the goal).

	let start_coords = L.A.pointToCoordinateArray(start_lat_lng),
	goal_coords = L.A.pointToCoordinateArray(goal_lat_lng),
	street_path_unordered = L.A.reversedCoordinates(lineSlice(start_coords, goal_coords, street_feature).geometry.coordinates);
	let start_to_path_beginning = start_lat_lng.distanceTo(L.latLng(street_path_unordered[0])),
	start_to_path_end = start_lat_lng.distanceTo(L.latLng(street_path_unordered[street_path_unordered.length - 1]));
	let street_path = start_to_path_beginning < start_to_path_end ?	street_path_unordered :	street_path_unordered.reverse();
	let street_path_lat_lngs = street_path.map(coords => {
		let lat_lng = L.latLng(coords);
		lat_lng.new_place = { type: "street", id: street_id },
		lat_lng.speed = speed;

		return lat_lng;
	});

	let first_lat = street_path_lat_lngs[0].lat,
	first_lng = street_path_lat_lngs[0].lng;

	//Exclude the last point if it's the same as the second to last point of this proposed segment,
	//and the second of it's the same as the first.
	//(since lineSlice adds a point for each other street in an intersection).
	if (street_path_lat_lngs.length > 1) {
		let second_lat = street_path_lat_lngs[1].lat,
		second_lng = street_path_lat_lngs[1].lng,
		final_lat = street_path_lat_lngs[street_path_lat_lngs.length - 1].lat,
		final_lng = street_path_lat_lngs[street_path_lat_lngs.length - 1].lng,
		penultimate_lat = street_path_lat_lngs[street_path_lat_lngs.length - 2].lat,
		penultimate_lng = street_path_lat_lngs[street_path_lat_lngs.length - 2].lng;

		if (first_lat === second_lat && first_lng === second_lng) {
			street_path_lat_lngs.shift();
		}

		if (final_lat === penultimate_lat && final_lng === penultimate_lng) {
			street_path_lat_lngs.pop();
		}
	}

	//Exclude the first point if it's already the last point of the already scheduled path.
	if (this.trip.path.length > 0) {
		let prev_lat = this.trip.path[this.trip.path.length - 1].lat,
		prev_lng = this.trip.path[this.trip.path.length - 1].lng;

		if (prev_lat === first_lat && prev_lng === first_lng) {
			street_path_lat_lngs.shift();
		}
	}

	this.trip.path.push(...street_path_lat_lngs);
}

/**
 * Schedule the agent up to travel between two points on a street network.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {LatLng} start_lat_lng - The coordinates of the point on the street from which the agent will be traveling.
 * @param {LatLng} goal_lat_lng - The coordinates of the point on the street to which the agent should travel.
 * @param {LatLng} start_int_lat_lng - The coordinates of the nearest intersection on the same street at the start_lat_lng.
 * @param {LatLng} goal_int_lat_lng - The coordinates of the nearest intersection on the same street as the goal_lat_lng.
 * @param {number} speed - The speed that the agent should travel, in meters per tick.
 */
Agent.setTravelOnStreetNetwork = function(start_lat_lng, goal_lat_lng, start_int_lat_lng, goal_int_lat_lng, speed) {
	let path = this.agentmap.getPath(start_int_lat_lng, goal_int_lat_lng, start_lat_lng, goal_lat_lng, true);

	for (let i = 0; i <= path.length - 2; i++) {
		let current_street_id = path[i].new_place.id,
		current_street_feature = this.agentmap.streets.getLayer(current_street_id).feature;

		this.setTravelOnSameStreet(path[i], path[i + 1], current_street_feature, current_street_id, speed);
	}
}

/**
 * Set a new, constant speed for the agent to move along its currently scheduled path.
 * @memberof Agent
 * @instance
 *
 * @param {number} speed - The speed (in meters per tick) that the agent should move. Must be >= .1.
 */
Agent.setSpeed = function(speed) {
	this.checkSpeed(speed);

	if (this.trip.goal_point !== null) {
		this.trip.speed = speed;
	}

	for (let spot of this.trip.path) {
		this.trip.speed = speed;
		spot.speed = speed;
	}
}

/**
 * Multiply the speed the agent moves along its currently scheduled path by a constant.
 * @memberof Agent
 * @instance
 *
 * @param {number} multiplier - The number to multiply the agent's scheduled speed by.
 * All scheduled speeds must be >= .1.
 */
Agent.multiplySpeed = function(multiplier) {
	if (this.trip.goal_point !== null) {
		this.trip.speed *= multiplier;
		this.checkSpeed(this.trip.speed);
	}

	for (let spot of this.trip.path) {
		spot.speed *= multiplier;
		this.checkSpeed(spot.speed);
	}
}

/**
 * Increase the speed the agent moves along its currently scheduled path by a constant.
 * @memberof Agent
 * @instance
 *
 * @param {number} magnitude - The number to add to the agent's scheduled speed.
 * All scheduled speeds must be >= .1
 */
Agent.increaseSpeed = function(magnitude) {
	if (this.trip.goal_point !== null) {
		this.trip.speed += magnitude;
		this.checkSpeed(this.trip.speed);
	}

	for (let spot of this.trip.path) {
		spot.speed += magnitude;
		this.checkSpeed(spot.speed);
	}
}

/**
 * Check whether a given speed is greater than the minimum.
 * @memberof Agent
 * @instance
 *
 * @param {number} speed - A number representing the speed of an agent in meters per second.
 */
Agent.checkSpeed = function(speed) {
	if (speed < .1) {
		throw new Error("Cannot assign speed below .1 to agent!");
	}
}

/**
 * Continue to move the agent directly along the points in its path, at approximately the speed associated with each point in the path.
 * Since two points along the path may be far apart, the agent will make multiple intermediary movements too, splitting up its transfer
 * from its current point to its goal point into a sub-path with multiple sub-goals.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {number} override_speed - Have the agent step this distance, instead of the distance suggested by the current state's speed property.
 */
Agent.travel = function(override_speed) {
	let current_coords = L.A.pointToCoordinateArray(this.trip.current_point),
	sub_goal_distance = override_speed ||this.trip.speed,
	sub_goal_coords = destination(current_coords, sub_goal_distance * .001,this.trip.angle).geometry.coordinates,
	sub_goal_lat_lng = L.latLng(L.A.reversedCoordinates(sub_goal_coords));

	let segment_to_goal = lineString([this.trip.current_point, this.trip.goal_point].map(point => L.A.pointToCoordinateArray(point))),
	segment_to_sub_goal = lineString([this.trip.current_point, sub_goal_lat_lng].map(point => L.A.pointToCoordinateArray(point)));

	let goal_lat_dist = Math.abs(this.trip.current_point.lat - this.trip.goal_point.lat),
	goal_lng_dist = Math.abs(this.trip.current_point.lng - this.trip.goal_point.lng);

	let dist_to_goal = length(segment_to_goal) * 1000,
	dist_to_sub_goal = length(segment_to_sub_goal) * 1000,
	leftover_after_goal;

	//Check if the distance to the sub_goal is greater than the distance to the goal, and if so, make the sub_goal equal the goal
	//and change the number of meters to the sub_goal to the number of meters to the goal.
	if (dist_to_goal < dist_to_sub_goal) {
		sub_goal_lat_lng = this.trip.goal_point,
		sub_goal_distance = dist_to_goal,
		leftover_after_goal = dist_to_sub_goal - dist_to_goal;
	}

	if (this.checkArrival(sub_goal_lat_lng, leftover_after_goal)) {
		return;
	}

	//Lat/Lng distance between current point and sub_goal point.
	let sub_goal_lat_dist = Math.abs(sub_goal_lat_lng.lat - this.trip.current_point.lat),
	sub_goal_lng_dist = Math.abs(sub_goal_lat_lng.lng - this.trip.current_point.lng);

	let half_meters = sub_goal_distance * 2,
	int_half_meters = Math.floor(half_meters),
	int_lat_step_value = this.trip.lat_dir * (sub_goal_lat_dist / half_meters),
	int_lng_step_value = this.trip.lng_dir * (sub_goal_lng_dist / half_meters),
	final_lat_step_value = this.trip.lat_dir * (sub_goal_lat_dist - Math.abs(int_lat_step_value * int_half_meters)),
	final_lng_step_value = this.trip.lng_dir * (sub_goal_lng_dist - Math.abs(int_lng_step_value * int_half_meters));

	//Intermediary movements.
	for (let i = 0; i < int_half_meters; ++i) {
		this.step(int_lat_step_value, int_lng_step_value);

		//If the agent is moving directly from a large distance, redirect it back towards the goal if it appears off course.
		if (this.trip.goal_point.move_directly === true) {
			let new_goal_lat_dist = Math.abs(this.trip.current_point.lat - this.trip.goal_point.lat),
			new_goal_lng_dist = Math.abs(this.trip.current_point.lng - this.trip.goal_point.lng);

			if (new_goal_lat_dist > goal_lat_dist || new_goal_lng_dist > goal_lng_dist) {
				this.travelTo(this.trip.goal_point);
			}
		}

		if (this.checkArrival(sub_goal_lat_lng, leftover_after_goal)) {
			return;
		}
	}

	//Last movement after intermediary movements.
	this.step(final_lat_step_value, final_lng_step_value, true);

	if (this.checkArrival(sub_goal_lat_lng, leftover_after_goal)) {
		return;
	}
};

/**
 * Move the agent a given latitude and longitude.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {number} lat_step_value - The number to add to the agent's latitude.
 * @param {number} lng_step_value - The number to add to the agent's longitude.
 */
Agent.step = function(lat_step_value, lng_step_value) {
	let new_lat_lng = L.latLng([this.trip.current_point.lat + lat_step_value, this.trip.current_point.lng + lng_step_value]);

	this.trip.current_point = new_lat_lng,
	this.steps_made++;

	//Only redraw the Agent's position if the number of steps the agent has moved is a multiple of the agentmap.animation_interval.
	if (this.agentmap.animation_interval > 0 && this.steps_made % this.agentmap.animation_interval === 0) {
		this.setLatLng(new_lat_lng);
	}
	else {
		this._latlng = new_lat_lng;
	}
};

/**
 * Check if the agent has arrived at the next goal in its path or to a sub_goal along the way and perform appropriate arrival operations.
 * @memberof Agent
 * @instance
 * @private
 *
 * @param {LatLng} sub_goal_lat_lng - A sub_goal on the way to the goal (possibly the goal itself).
 * @param {number} leftover_after_goal - If the agent arrives at its goal during the tick, the number of meters, according to its speed,
 * leftover beyond the goal that it should still move during the tick.
 */
Agent.checkArrival = function(sub_goal_lat_lng, leftover_after_goal) {
	if (this.trip.goal_point.distanceTo(this.trip.current_point) < .1) {
		this.place = this.trip.path[0].new_place;
		arrived = true;

		this.trip.path.shift();

		if (this.trip.path.length === 0) {
			this.resetTrip();
		}
		else {
			this.travelTo(this.trip.path[0]);

			//If it still needs to move a certain distance during this tick, move it that distance towards the next goal before returning.
			if (leftover_after_goal > 0) {
				this.travel(leftover_after_goal);
			}
		}

		this.trip.moving = false;

		return true;
	}
	else if (sub_goal_lat_lng.distanceTo(this.trip.current_point) < .1) {
		this.trip.moving = false;

		return true;
	}
};

/**
 * Make the agent proceed along its trip.
 * @memberof Agent
 * @instance
 */
Agent.moveIt = function() {
	//Make sure the agent isn't paused or already moving.
	if (!this.trip.paused && !this.trip.moving) {
		//Call the agent's fine_controller before it begins moving.
		this.fine_controller();

		//Check if the agent has a goal point, and if so travel towards it.
		if (this.trip.goal_point !== null) {
			this.trip.moving = true;
			this.travel();
		}
		//Otherwise, if there's a scheduled path that the agent hasn't started traveling on yet,
		//start traveling on it.
		else if (this.trip.path.length !== 0) {
			this.trip.moving = true;
			this.startTrip();
			this.travel();
		}
	}
}

Agent = L.CircleMarker.extend(Agent);

/**
 * Returns an agent object.
 *
 * @param {LatLng} lat_lng - A pair of coordinates to locate the agent at.
 * @param {Object} options - An array of options for the agent, namely its layer.
 * @param {Agentmap} agentmap - The agentmap instance in which the agent exists.
 */
function agent(lat_lng, options, agentmap) {
	return new Agent(lat_lng, options, agentmap);
}

/**
 * A user-defined callback function that returns a feature with appropriate geometry and properties to represent an agent.
 *
 * @callback agentFeatureMaker
 * @param {number} id - The agent's Leaflet layer ID.
 * @returns {Point} - a GeoJSON Point feature with properties and coordinates for the agent, including
 * a "place" property that will set the agent's initial {@link Place} and an object "layer_options" property
 * that will specify the feature's Leaflet options (like its color, size, etc.). All other provided properties
 * will be transferred to the Agent object once it is created.
 * See {@link https://leafletjs.com/reference-1.3.2.html#circlemarker} for all possible layer options.
 *
 * @example
 * let point = {
 * 	"type": "Feature",
 * 	"properties": {
 * 		"layer_options": {
 * 			"color": "red",
 * 			"radius": .5,
 * 		},
 * 		"place": {
 * 			"type": "unit",
 * 			"id": 89
 * 		},
 *
 * 		age: 72,
 * 		home_city: "LA"
 * 	},
 * 	"geometry" {
 * 		"type": "Point",
 * 		"coordinates": [
 * 			14.54589,
 * 			57.136239
 * 		]
 * 	}
 * }
 */

/**
 * A standard {@link agentFeatureMaker}, which sets an agent's location to be the point near the center of the iᵗʰ unit of the map,
 * its place property to be that unit's, and its layer_options to be red and of radius .5 meters.
 * @memberof Agentmap
 * @instance
 * @type {agentFeatureMaker}
 */
function seqUnitAgentMaker(id){
	let index = this.agents.count();

	if (index > this.units.getLayers().length - 1) {
		throw new Error("seqUnitAgentMaker cannot accommodate more agents than there are units.");
	}

	let unit = this.units.getLayers()[index],
	unit_id = this.units.getLayerId(unit),
	center_point = centroid(unit.feature);
	center_point.properties.place = {"type": "unit", "id": unit_id},
	center_point.properties.layer_options = {radius: .5, color: "red", fillColor: "red"};

	return center_point;
}

/**
 * Generate some number of agents and place them on the map.
 * @memberof Agentmap
 * @instance
 *
 * @param {number} count - The desired number of agents.
 * @param {agentFeatureMaker} agentFeatureMaker - A callback that determines an agent i's feature properties and geometry (always a Point).
 */
function agentify(count, agentFeatureMaker) {
	let agentmap = this;

	if (!(this.agents instanceof L.LayerGroup)) {
		this.agents = L.featureGroup().addTo(this.map);
	}

	let agents_existing = agentmap.agents.getLayers().length;
	for (let i = agents_existing; i < agents_existing + count; i++) {
		let new_agent = agent(null, null, agentmap);

		//Callback function aren't automatically bound to the agentmap.
		let boundFeatureMaker = agentFeatureMaker.bind(agentmap),
		agent_feature = boundFeatureMaker(new_agent._leaflet_id);

		let coordinates = L.A.reversedCoordinates(agent_feature.geometry.coordinates),
		place = agent_feature.properties.place,
		layer_options = agent_feature.properties.layer_options;

		//Make sure the agent feature is valid and has everything we need.
		if (!L.A.isPointCoordinates(coordinates)) {
			throw new Error("Invalid feature returned from agentFeatureMaker: geometry.coordinates must be a 2-element array of numbers.");
		}
		else if (typeof(place.id) !== "number") {
			throw new Error("Invalid feature returned from agentFeatureMaker: properties.place must be a {unit: unit_id} or {street: street_id} with an existing layer's ID.");
		}

		new_agent.setLatLng(coordinates);
		new_agent.setStyle(layer_options);

		delete agent_feature.properties.layer_options;
		Object.assign(new_agent, agent_feature.properties);

		this.agents.addLayer(new_agent);
	}
}

Agentmap.prototype.agent = agent,
Agentmap.prototype.agentify = agentify,
Agentmap.prototype.seqUnitAgentMaker = seqUnitAgentMaker;

exports.Agent = Agent,
exports.agent = agent;

},{"./agentmap":43,"./routing":47,"@turf/along":2,"@turf/bearing":4,"@turf/boolean-point-in-polygon":5,"@turf/buffer":6,"@turf/centroid":10,"@turf/destination":13,"@turf/helpers":15,"@turf/length":18,"@turf/line-slice":24,"@turf/nearest-point-on-line":36}],45:[function(require,module,exports){
/* Functions that help design and generate building units onto the map. */

let bearing = require('@turf/bearing').default,
destination = require('@turf/destination').default,
along = require('@turf/along').default,
lineIntersect = require('@turf/line-intersect').default,
intersect = require('@turf/intersect').default,
Agentmap = require('./agentmap').Agentmap,
streetsToGraph = require('./routing').streetsToGraph,
getPathFinder = require('./routing').getPathFinder;

/**
 * Generate and setup the desired map features (e.g. streets, houses).
 * @memberof Agentmap
 * @instance
 *
 * @param {Array.<Array.<number>>} bounding_box - The map's top-left and bottom-right coordinates.
 * @param {object} streets_data - A GeoJSON Feature Collection object containing the OSM street features inside the bounding box.
 * @param {object} [street_options] - An object containing the Leaflet styling options for streets. See available options here: {@link https://leafletjs.com/reference-1.3.2.html#polyline-l-polyline}.
 * @param {object} [unit_options] - An object containing the Leaflet & AgentMaps styling options for units.<br/>See available Leaflet options here: {@link https://leafletjs.com/reference-1.3.2.html#polygon-l-polygon}<br/>Additional AgentMaps-specific options are described below.
 * @param {number} [unit_options.front_buffer = 6] - The number of meters beetween the front of unit and its street.
 * @param {number} [unit_options.side_buffer = 3] - The number of meters between two units on the same street.
 * @param {number} [unit_options.length = 14] - The length of the unit in meters along the street.
 * @param {number} [unit_options.depth = 18] - The depth of the unit in meters out from its front.
 * @param {object} [units_data]- If you want to load a previously generated AgentMaps.units object instead of generating one from scarch: A GeoJSON Feature Collection of an AgentMaps.units featureGroup.
 */
function buildingify(bounding_box, streets_data, street_options, unit_options, units_data) {
	setupStreetFeatures.call(this, streets_data, street_options);
	setupUnitFeatures.call(this, bounding_box, streets_data, unit_options, units_data);
}

/**
 * Generate and setup streets based on the provided GeoJSON data.
 *
 * @param {object} streets_data - A GeoJSON Feature Collection object containing the OSM street features inside the bounding box.
 * @param {object} street_options - An object containing the Leaflet styling options for streets.
 */
function setupStreetFeatures(streets_data, street_options) {
	let street_features = getStreetFeatures(streets_data);
	
	default_options = {
		"color": "yellow",
		"weight": 4,
		"opacity": .5
	};

	street_options = Object.assign(default_options, street_options);

	let street_feature_collection = {
		type: "FeatureCollection",
		features: street_features
	};

	this.streets = L.geoJSON(
		street_feature_collection,
		street_options
	).addTo(this.map);

	//Map streets' OSM IDs to their Leaflet IDs.
	this.streets.id_map = {};
	
	//Having added the streets as layers to the map, do any processing that requires access to those layers.
	this.streets.eachLayer(function(street) {
		this.streets.id_map[street.feature.id] = street._leaflet_id; 
		
		addStreetLayerIntersections.call(this, street);
	}, this);
	
	this.streets.graph = streetsToGraph(this.streets),
	this.pathfinder = getPathFinder(this.streets.graph);
}

/**
 * Get all streets from the GeoJSON data.
 * @private
 *
 * @param {Object} streets_data - A GeoJSON Feature Collection object containing the OSM streets inside the bounding box.
 * @returns {Array<Feature>} -  array of street features.
 */
function getStreetFeatures(streets_data) {
	let street_features = [];

	for (let i =  0; i < streets_data.features.length; ++i) {
		let feature = streets_data.features[i];
		
		if (feature.geometry.type === "LineString" && feature.properties.tags.highway) {
			let street_feature = feature;

			street_features.push(street_feature);
		}
	}

	return street_features;
}

/**
 * Gets the intersections of all the streets on the map and adds them as properties to the street layers.
 * @private
 * 
 * @param {object} street - A Leaflet polyline representing a street.
 */
function addStreetLayerIntersections(street) {
	let street_id = street._leaflet_id;

	street.intersections = typeof(street.intersections) === "undefined" ? {} : street.intersections;

	this.streets.eachLayer(function(other_street) {
		let other_street_id = other_street._leaflet_id;

		//Skip if both streets are the same, or if the street already has its intersections with the other street.
		if (typeof(street.intersections[other_street_id]) === "undefined" && street_id !== other_street_id) {
			let street_coords = street.getLatLngs().map(L.A.pointToCoordinateArray),
			other_street_coords = other_street.getLatLngs().map(L.A.pointToCoordinateArray),
			identified_intersections = L.A.getIntersections(street_coords, other_street_coords, [street_id, other_street_id]).map(
				identified_intersection => 
				[L.latLng(L.A.reversedCoordinates(identified_intersection[0])), identified_intersection[1]]
			);

			if (identified_intersections.length > 0) {
				street.intersections[other_street_id] = identified_intersections,
				other_street.intersections = typeof(other_street.intersections) === "undefined" ? {} : other_street.intersections,
				other_street.intersections[street_id] = identified_intersections;
			}
		}
	});
}

/**
 * Generate and setup building units based on the provided GeoJSON data.
 *
 * @param {Array.<Array.<number>>} bounding_box - The map's top-left and bottom-right coordinates.
 * @param {object} streets_data - A GeoJSON Feature Collection object containing the OSM street features inside the bounding box.
 * @param {object} unit_options - An object containing the Leaflet & AgentMaps styling options for units.
 * @param {object} [units_data] - If you want to load a previously generated AgentMaps.units object instead of generating one from scarch: A GeoJSON Feature Collection of an AgentMaps.units featureGroup.
 */
function setupUnitFeatures(bounding_box, streets_data, unit_options = {}, units_data) {
	let default_options = {
			"color": "green",
			"weight": 1,
			"opacity": .87,
			"front_buffer": 6,
			"side_buffer": 3,
			"length": 14,
			"depth": 18
	};

	unit_options = Object.assign(default_options, unit_options);
	
	let unit_feature_collection;

	//If no units_data is supplied, generate the units from scratch.
	if (typeof(units_data) === "undefined") {
		//Bind getUnitFeatures to "this" so it can access the agentmap as "this.agentmap".
		let unit_features = getUnitFeatures.bind(this)(bounding_box, streets_data, unit_options);

		unit_feature_collection = { 
			type: "FeatureCollection", 
			features: unit_features
		};
	}
	else {
		unit_feature_collection = units_data;
	}
	
	this.units = L.geoJSON(
		unit_feature_collection,
		unit_options
	).addTo(this.map);

	//Having added the units as layers to the map, do any processing that requires access to those layers.
	this.units.eachLayer(function(unit) {
		if (typeof(units_data) === "undefined") {
			unit.street_id = unit.feature.properties.street_id;
		}
		else {
			unit.street_id = this.streets.id_map[unit.feature.properties.OSM_street_id];
		}

		unit.street_anchors = unit.feature.properties.street_anchors,
		//Change the IDs of each unit in this unit's neighbours array into the appropriate Leaflet IDs.
		unit.neighbors = getUnitNeighborLayerIDs.call(this, unit.feature.properties.neighbors);
	}, this);
}

/**
 * Given an array of pre-layer IDs, check if any of them correspond to the pre-layer IDs of unit layers, and if so
 * return an array of the corresponding layer IDs.
 * @private
 *
 * @param {Array<?number>} - An array of pre-layer feature IDs for a unit's neighbors.
 * @returns {Array<?number>} - An array of Leaflet layer IDs corresponding to the unit's neighbors.
 */
function getUnitNeighborLayerIDs(neighbors) {
	let neighbor_layer_ids = neighbors.map(function(neighbor) {
		if (neighbor !== null) {
			let neighbor_layer_id = null;
			
			this.units.eachLayer(function(possible_neighbor_layer) {
				if (possible_neighbor_layer.feature.properties.id === neighbor) {
					neighbor_layer_id = this.units.getLayerId(possible_neighbor_layer);
				}
			}, this);

			return neighbor_layer_id;
		}
		else {
			return null;
		}
	}, this);

	return neighbor_layer_ids;
}

/**
 * Get all appropriate units within the desired bounding box.
 * @private
 *
 * @param {Array.<Array.<number>>} bounding_box - The map's top-left and bottom-right coordinates.
 * @param {Object} streets_data - A GeoJSON Feature Collection object containing the OSM street features inside the bounding box.
 * @param {object} unit_options - An object containing the AgentMaps styling options for units.
 * @returns {Array<Feature>} -  array of features representing real estate units.
 */
function getUnitFeatures(bounding_box, streets_data, unit_options) {
	let proposed_unit_features = [];
	
	this.streets.eachLayer(function(layer) {
		let street_feature = layer.feature,
		street_id = layer._leaflet_id,
		street_OSM_id = layer.feature.id,
		proposed_anchors = getUnitAnchors(street_feature, bounding_box, unit_options),
		new_proposed_unit_features = generateUnitFeatures(proposed_anchors, proposed_unit_features, street_id, street_OSM_id, unit_options);
		proposed_unit_features.push(...new_proposed_unit_features);
	});

	unit_features = unitsOutOfStreets(proposed_unit_features, this.streets);
	
	return unit_features;
}

/**
 * Given an array of anchor pairs, for each anchor pair find four 
 * nearby points on either side of the street appropriate to build a unit(s) on.
 * @private
 *
 * @param {Array<Array<Feature>>} unit_anchors - Array of pairs of points around which to anchor units along a street.
 * @param {Array<Feature>} proposed_unit_features - Array of features representing building units already proposed for construction.
 * @param {string} street_leaflet_id - The Leaflet layer ID of the street feature along which the unit is being constructed.
 * @param {string} street_OSM_id - The OSM feature ID of the street feature along which the unit is being constructed.
 * @param {object} unit_options - An object containing the AgentMaps styling options for units.
 * @returns {Array<Feature>} unit_features - Array of features representing units.
 */
function generateUnitFeatures(unit_anchors, proposed_unit_features, street_leaflet_id, street_OSM_id, unit_options) {
	//One sub-array of unit features for each side of the road.
	let unit_features = [[],[]],
	starting_id = proposed_unit_features.length,
	increment = 1;
	
	for (let anchor_pair of unit_anchors) {
		//Pair of unit_features opposite each other on a street.
		let unit_pair = [null, null];
		
		for (let i of [1, -1]) {
			let anchor_a = anchor_pair[0].geometry.coordinates,
			anchor_b = anchor_pair[1].geometry.coordinates,
			anchor_latLng_pair = [anchor_a, anchor_b],
			street_buffer = unit_options.front_buffer / 1000, //Distance between center of street and start of unit.
			house_depth = unit_options.depth / 1000,
			angle = bearing(anchor_a, anchor_b),
			new_angle = angle + i * 90, //Angle of line perpendicular to the anchor segment.
			unit_feature = { 
				type: "Feature",
				properties: {
					street: "none"
				},
				geometry: {
					type: "Polygon",
					coordinates: [[]]
				}
			};
			unit_feature.geometry.coordinates[0][0] = destination(anchor_a, street_buffer, new_angle).geometry.coordinates,
			unit_feature.geometry.coordinates[0][1] = destination(anchor_b, street_buffer, new_angle).geometry.coordinates,
			unit_feature.geometry.coordinates[0][2] = destination(anchor_b, street_buffer + house_depth, new_angle).geometry.coordinates,
			unit_feature.geometry.coordinates[0][3] = destination(anchor_a, street_buffer + house_depth, new_angle).geometry.coordinates;
			unit_feature.geometry.coordinates[0][4] = unit_feature.geometry.coordinates[0][0];

			//Exclude the unit if it overlaps with any of the other proposed units.
			let all_proposed_unit_features = unit_features[0].concat(unit_features[1]).concat(proposed_unit_features);
			if (noOverlaps(unit_feature, all_proposed_unit_features)) { 
				//Recode index so that it's useful here.
				i = i === 1 ? 0 : 1;

				unit_feature.properties.street_id = street_leaflet_id,
				unit_feature.properties.OSM_street_id = street_OSM_id,
				unit_feature.properties.street_anchors = anchor_latLng_pair,	
				unit_feature.properties.neighbors = [null, null, null],
				unit_feature.properties.id = starting_id + increment,
				increment += 1;
				
				if (unit_features[i].length !== 0) {
					//Make previous unit_feature this unit_feature's first neighbor.
					unit_feature.properties.neighbors[0] = unit_features[i][unit_features[i].length - 1].properties.id,
					//Make this unit_feature the previous unit_feature's second neighbor.
					unit_features[i][unit_features[i].length - 1].properties.neighbors[1] = unit_feature.properties.id;
				}
				
				if (i === 0) {
					unit_pair[0] = unit_feature;
				}
				else {
					if (unit_pair[0] !== null) {
						//Make unit_feature opposite to this unit_feature on the street its third neighbor.
						unit_feature.properties.neighbors[2] = unit_pair[0].properties.id,
						//Make unit_feature opposite to this unit_feature on the street's third neighbor this unit_feature.
						unit_pair[0].properties.neighbors[2] = unit_feature.properties.id;
					}
					
					unit_pair[1] = unit_feature;
				}
			}
		}
		
		if (unit_pair[0] !== null) {
			unit_features[0].push(unit_pair[0]);
		}

		if (unit_pair[1] !== null) {
			unit_features[1].push(unit_pair[1]);
		}
	}

	let unit_features_merged = [].concat(...unit_features);

	return unit_features_merged;
}

/**
 * Find anchors for potential units. chors are the pairs of start 
 * and end points along the street from which units may be constructed.
 * @private
 * 
 * @param {Feature} street_feature - A GeoJSON feature object representing a street.
 * @param {object} unit_options - An object containing the AgentMaps styling options for units.
 * @returns {Array<Array<Feature>>} - Array of pairs of points around which to anchor units along a street.  
 */
function getUnitAnchors(street_feature, bounding_box, unit_options) {
	let unit_anchors = [],
	unit_length = unit_options.length / 1000, //Kilometers.
	unit_buffer = unit_options.side_buffer / 1000, //Distance between units, kilometers.
	endpoint = street_feature.geometry.coordinates[street_feature.geometry.coordinates.length - 1],
	start_anchor = along(street_feature, 0),
	end_anchor = along(street_feature, unit_length),
	distance_along = unit_length;

	while (end_anchor.geometry.coordinates != endpoint) {
		//Exclude proposed anchors if they're outside of the bounding box.
		start_coord = L.A.reversedCoordinates(start_anchor.geometry.coordinates), 
		end_coord = L.A.reversedCoordinates(end_anchor.geometry.coordinates);
		if (L.latLngBounds(bounding_box).contains(start_coord) &&
			L.latLngBounds(bounding_box).contains(end_coord)) {
				unit_anchors.push([start_anchor, end_anchor]);
		}

		//Find next pair of anchors.
		start_anchor = along(street_feature, distance_along + unit_buffer);
		end_anchor = along(street_feature, distance_along + unit_buffer + unit_length);
		
		distance_along += unit_buffer + unit_length
	}

	return unit_anchors;
}

/**
 * Get an array of units excluding units that overlap with streets.
 * @private
 *
 * @param {Array<Feature>} unit_features - Array of features representing units.
 * @param {Array<Layer>} street_layers - Array of Leaflet layers representing streets.
 * @returns {Array<Feature>} - unit_features, but with all units that intersect any streets removed.
 */
function unitsOutOfStreets(unit_features, street_layers) {
	let processed_unit_features = unit_features.slice();
	
	street_layers.eachLayer(function(street_layer) {
		let street_feature = street_layer.feature;
		for (let unit_feature of processed_unit_features) {
			let intersection_exists = lineIntersect(street_feature, unit_feature).features.length > 0;
			if (intersection_exists) {
				processed_unit_features.splice(processed_unit_features.indexOf(unit_feature), 1, null);
			}
		}	
	
		processed_unit_features = processed_unit_features.filter(feature => feature === null ? false : true);
	});
	

	return processed_unit_features;
}

/**
 * Check whether a polygon overlaps with any member of an array of polygons.
 * @private
 *
 * @param {Feature} reference_polygon_feature - A geoJSON polygon feature.
 * @param {Array<Feature>} polygon_feature_array - Array of geoJSON polygon features.
 * @returns {boolean} - Whether the polygon_feature overlaps with any one in the array.
 */	
function noOverlaps(reference_polygon_feature, polygon_feature_array) {
	//return true;
	for (feature_array_element of polygon_feature_array) {
		let overlap_exists = intersect(reference_polygon_feature, feature_array_element);
		if (overlap_exists) {
			return false;
		}
	}

	return true;
}

Agentmap.prototype.buildingify = buildingify;

},{"./agentmap":43,"./routing":47,"@turf/along":2,"@turf/bearing":4,"@turf/destination":13,"@turf/intersect":16,"@turf/line-intersect":20}],46:[function(require,module,exports){
let agentmap = require('./agentmap'),
agents = require('./agents'),
buildings = require('./buildings'),
utils = require('./utils');

L.A = Object.assign({}, agentmap, agents, buildings, utils);

},{"./agentmap":43,"./agents":44,"./buildings":45,"./utils":48}],47:[function(require,module,exports){
//** Here we have utilities to convert OSM geojson data into a distance-weighted graph and find the shortest path between two points. **//

let path = require("ngraph.path"),
createGraph = require("ngraph.graph"),
lineSlice = require('@turf/line-slice').default,
length = require('@turf/length').default,
Agentmap = require('./agentmap').Agentmap;

/**
 * Convert a layerGroup of streets into a graph.
 * @private
 *
 * @param {LayerGroup} streets - A Leaflet layerGroup of streets, forming a street network.
 * @returns {Object} - A graph representing the street network, operable by the ngraph pathfinder. 
 */
function streetsToGraph(streets) {
	let graph = createGraph(),
	streetToGraphBound = streetToGraph.bind(this, graph);
	
	//For each street, get an array of indices for the start, intersections, and end coordinates, in order from
	//start to end. Then, add the coordinates at each index as a node, and an edge between each adjacent node in the array,
	//associating the distance between the nodes (between their coordinates) with each edge.
	streets.eachLayer(streetToGraphBound);

	return graph;
}

/**
 * Process a street layer and add it into a graph.
 *
 */
function streetToGraph(graph, street) {
	let street_id = street._leaflet_id,
	intersection_indices = [],
	street_points = street.getLatLngs();
	
	//Populate intersection_indices with the indices of all of the street's intersections in its coordinate array.
	for (let cross_street in street.intersections) {
		let intersections = street.intersections[cross_street];
		
		for (let intersection of intersections) {
			let intersection_index = intersection[1][street_id];
			
			//Ignore duplicate intersection points (caused by 3-way intersections).
			if (!intersection_indices.some(other_intersection_index => other_intersection_index === intersection_index)) {
				intersection_indices.push(intersection_index);
			}
		}
	}

	//Sort the intersection_indices so that they are in order from the start of the street's coordinate array to the end;
	//this is why we're not getting the raw coordinates, but their indices first, so they can be sorted.
	intersection_indices = intersection_indices.sort(function(a, b) {
		return a - b;
	});

	//Check if beginning and end points of the street are in the intersection_incides; if not, add them.
	if (!intersection_indices.some(intersection_index => intersection_index === 0)) {
		intersection_indices.unshift(0);
	}
	if (!intersection_indices.some(intersection_index => intersection_index === street_points.length - 1)) {
		intersection_indices.push(street_points.length - 1);
	}

	//Make a graph out of segments of the street between the start, intersections, and end of the street,
	//so that the nodes are the coordinates of the start, end, and intersection points, and the edges are
	//the segments between successive nodes. Each edge is associated with the geographic distance between its nodes.
	for (let i = 0; i <= intersection_indices.length - 2; i++) {
		let node_a = street_points[intersection_indices[i]],
		node_b = street_points[intersection_indices[i + 1]],
		a_string = encodeLatLng(node_a),
		b_string = encodeLatLng(node_b),
		start_coords = L.A.pointToCoordinateArray(node_a),
		end_coords = L.A.pointToCoordinateArray(node_b),
		segment = lineSlice(start_coords, end_coords, street.toGeoJSON()),
		distance = length(segment);
		graph.addLink(a_string, b_string, {
			distance: distance,
			place: { type: "street",
				id: street_id } 
		});
	}
}

/**
 * Given an OSM street network (graph), return an A* pathfinder that can operate on it.
 * @private
 * 
 * @param {object} graph - An ngraph graph representing an OSM street network.
 * @returns {object} - An A* pathfinder for the graph.
 */
function getPathFinder(graph) {
	return path.aStar(graph, {
		distance(fromNode, toNode, link) {
			return link.data.distance;
		}
	});
}

/**
 * Get a path between two points on a graph.
 * @memberof Agentmap
 * @instance
 * @private
 *
 * @param start_int_lat_lng {LatLng} - The coordinates of the nearest intersection on the same street at the start_lat_lng.
 * @param goal_int_lat_lng {LatLng} - The coordinates of the nearest intersection on the same street as the goal_lat_lng.
 * @param start_lat_lng {LatLng} - The coordinates of the point on the street from which the agent will be traveling.
 * @param goal_lat_lng {LatLng} - The coordinates of the point on the street to which the agent should travel.
 * @param {Boolean} [sparse=false] - Whether to exclude intersections between the first and last along a street-specific path (which are superfluous for extracting the necessary sub-street).
 * @return {Array<Array<number>>} - An array of points along the graph, leading from the start to the end.
 */
function getPath(start_int_lat_lng, goal_int_lat_lng, start_lat_lng, goal_lat_lng, sparse = false) {
	let start_coord = encodeLatLng(start_int_lat_lng),
	end_coord = encodeLatLng(goal_int_lat_lng),
	encoded_path = this.pathfinder.find(start_coord, end_coord),
	path = [];
	
	if (encoded_path.length > 0 && decodeCoordString(encoded_path[0].id).distanceTo(start_int_lat_lng) > 
					decodeCoordString(encoded_path[0].id).distanceTo(goal_int_lat_lng)) {
		encoded_path = encoded_path.reverse();
	}

	if (sparse === true && encoded_path.length >= 2) {
		let sparse_path = [], 
		recent_street = null,
		current_street = null;
		
		for (let i = 0; i <= encoded_path.length - 2; i++) {
			current_street = this.streets.graph.getLink(encoded_path[i].id, encoded_path[i + 1].id) ||
				this.streets.graph.getLink(encoded_path[i + 1].id, encoded_path[i].id);
			
			if (recent_street === null || current_street.data.place.id !== recent_street.data.place.id) {
				let decoded_coords = decodeCoordString(encoded_path[i].id, current_street.data.place);
				sparse_path.push(decoded_coords);
			}
				
			//If the last place on the path to the goal is labeled with a different street id than the goal,
			//add it to the sparse path.	
			if (i === encoded_path.length - 2) {
				let decoded_coords = decodeCoordString(encoded_path[i + 1].id, current_street.data.place);
				sparse_path.push(decoded_coords);
			}
		}
			
		path = sparse_path;
	}
	else {
		path = encoded_path.map(point => decodeCoordString(point.id, 0));
	}
	
	path.unshift(start_lat_lng);
	path.push(goal_lat_lng);
	
	//If the goal point lies before the first intersection of the goal street, then the 2nd to last point in the
	//path will have the previous street's id attached to it. If the goal lies on a different street, make
	//sure the 2nd to last point (the street path intersection point before the goal) has the same street id as the goal.
	if (path[path.length - 2].new_place.id !== goal_lat_lng.new_place.id) {
		path[path.length - 2].new_place = goal_lat_lng.new_place;
	}

	//If the second [to last] point--namely the intersection closest to the start [goal]--is further from the third
	//[to last] point than the goal, and all three points are on the same street, remove the second [to last] point.
	if (path.length >= 3) {
		checkStartExcess.call(this, path);
		checkEndExcess.call(this, path);
	}
	
	return path;
}

//checkStartExcess and checkEndExcess are _much_ easier to follow given distinct variable names,
//and so they are not abstracted into one more general function.

/** 
 * If the first two points after the start point share the same street as the start point, and the
 * third point is closer to the first (start) point than it is to the second point, remove the 
 * second point, as it's a superfluous detour.<br/><br/>
 *
 * Typically happens when the start point's nearest intersection is beyond it on the street,
 * and so the path would have an agent travel from the start, then to the intersection,
 * then backwards to the third point.
 * @private
 *
 * @param {Array<LatLng>} path - An array of LatLngs representing a path for an agent to travel along.
 */
function checkStartExcess(path) {
	let start_street = this.streets.getLayer(path[0].new_place.id),
	second_street_id = path[1].new_place.id,	
	start_second_intersections = start_street.intersections[second_street_id],
	second_is_intersection = typeof(start_second_intersections) === "undefined" ? false :
		start_second_intersections.some(intersection => 
		intersection[0].lat === path[1].lat && intersection[0].lng === path[1].lng),
	third_street_id = path[2].new_place.id,
	start_third_intersections = start_street.intersections[third_street_id],
	third_is_intersection = typeof(start_third_intersections) === "undefined" ? false :
		start_third_intersections.some(intersection =>
		intersection[0].lat === path[2].lat && intersection[0].lng === path[2].lng);

	if ((second_is_intersection || second_street_id === path[0].new_place.id) && 
		(third_is_intersection || third_street_id === path[0].new_place.id)) {
		if (path[2].distanceTo(path[0]) <
			path[2].distanceTo(path[1])) {
			path.splice(1, 1);
		}
	}
}

/** 
 * If the last two points before the goal point share the same street as the goal point, and the
 * first point is closer to the third (goal) point than it is to the second point, remove the 
 * second point, as it's a superfluous detour.<br/><br/>
 *
 * Typically happens when the goal point's nearest intersection is beyond it on the street,
 * and so the path would have an agent travel from the first point, then to the intersection (second point),
 * then backwards to the (third) goal point.<br/><br/>
 *
 * @private
 *
 * @param {Array<LatLng>} path - An array of LatLngs representing a path for an agent to travel along.
 */
function checkEndExcess(path) {
	let goal_street = this.streets.getLayer(path[path.length - 1].new_place.id),
	second_to_last_street_id = path[path.length - 2].new_place.id,
	goal_second_to_last_intersections = goal_street.intersections[second_to_last_street_id],
	second_to_last_is_intersection = typeof(goal_second_to_last_intersections) === "undefined" ? false :
		goal_second_to_last_intersections.some(intersection => 
		intersection[0].lat === path[path.length - 1].lat && intersection[0].lng === path[path.length - 1].lng),
	third_last_street_id = path[path.length - 3].new_place.id,
	goal_third_last_intersections = goal_street.intersections[third_last_street_id],
	third_last_is_intersection = typeof(goal_third_last_intersections) === "undefined" ? false :
		goal_third_last_intersections.some(intersection =>
		intersection[0].lat === path[path.length - 3].lat && intersection[0].lng === path[path.length - 3].lng);

	if ((second_to_last_is_intersection || second_to_last_street_id === path[path.length - 1].new_place.id) &&
		(third_last_is_intersection || third_last_street_id === path[path.length - 1].new_place.id) && 
		path.length >= 3) {
		if (path[path.length - 3].distanceTo(path[path.length - 1]) <
			path[path.length - 3].distanceTo(path[path.length - 2])) {
			path.splice(path.length - 2, 1);
		}
	}
}

/**
 * Turn a LatLng object into a string representing its coordinates (to act as a graph node's ID).
 * @private
 *
 * @param {LatLng} lat_lng - The coordinates to encode into a string.
 * @returns {string} - A string containing coordinates in the format of "Latitude,Longitude".
 */
function encodeLatLng(lat_lng) {
	return lat_lng.lat.toString() + "," + lat_lng.lng.toString();
}

/**
 * Turn a string containing coordinates (a graph node's ID) into a LatLng object.
 * @private
 *
 * @param {string} coord_string - A string containing coordinates in the format of "Latitude,Longitude".
 * @param {object} place - An object specifying the place of the coordinate string.
 * @returns {LatLng} - The coordinates encoded by the coord_string.
 */
function decodeCoordString(coord_string, place) {
	let coord_strings = coord_string.split(","),
	lat_lng = L.latLng(coord_strings);
	lat_lng.new_place = place;

	return lat_lng;
}

Agentmap.prototype.getPath = getPath;

exports.streetsToGraph = streetsToGraph;
exports.getPathFinder = getPathFinder;
exports.encodeLatLng = encodeLatLng;

},{"./agentmap":43,"@turf/length":18,"@turf/line-slice":24,"ngraph.graph":55,"ngraph.path":64}],48:[function(require,module,exports){
/* A few functions that may be useful in other modules. */

/**
 * Given a geoJSON geometry object's coordinates, return the object, but with
 * all the coordinates reversed. <br /point.geometry && point.geometry.coordinates && >
 * 
 * Why? GeoJSON coordinates are in lngLat format by default, while Leaflet uses latLng.
 * L.geoJSON will auto-reverse the order of a GeoJSON object's coordinates, as it
 * expects geoJSON coordinates to be lngLat. However, normal, non-GeoJSON-specific Leaflet
 * methods expect Leaflet's latLng pairs and won't auto-reverse, so we have to do that
 * manually if we're preprocessing the GeoJSON data before passing it to L.geoJSON.
 * 
 * @param {Array<number|Array<number|Array<number>>>} coordinates - GeoJSON coordinates for a point, (multi-)line, or (multi-)polygon.
 * @returns {Array<number|Array<number|Array<number>>>} - Reversed geoJSON coordinates for a point, (multi-)line, or (multi-)polygon.
 */
function reversedCoordinates(coordinates) {
	let reversed = coordinates.slice();
	if (typeof coordinates[0] != "number") {
		for (let inner_coordinates of coordinates) {
			reversed.splice(reversed.indexOf(inner_coordinates), 1, reversedCoordinates(inner_coordinates));
		}
	}
	else {
		reversed = [coordinates[1], coordinates[0]];
	}

	return reversed;
}

/**
 * Given an array, check whether it can represent the coordinates of a point.
 *
 * @param {Array} array - Array to check.
 * @returns {boolean} - Whether the array can be the coordinates of a point.
 */
function isPointCoordinates(array) {
	if (array.length !== 2 || 
		typeof(array[0]) !== "number" ||
		typeof(array[1]) !== "number") {
		return false;
	}

	return true;
}

/**
 * Given either a GeoJSON feature, L.latLng, or coordinate array containing the coordinates of a point,
 * return an array of the coordinates.
 *
 * @param {Point|Array<number>|LatLng} point - The data containing the point's coordinates (latitude & longitude).
 * @returns {Array<number>} - Array of the point's coordinates. I.e.: [lng, lat].
 */
function pointToCoordinateArray(point) {
	let coordinate_array;

	if (typeof(point.lat) === "number" && typeof(point.lng) === "number") {
		coordinate_array = [point.lng, point.lat];
	}
	else if (point.geometry && point.geometry.coordinates && isPointCoordinates(point.geometry.coordinates)) {
		coordinate_array = point.geometry.coordinates;
	}
	else if (isPointCoordinates(point)) {
		coordinate_array = point;
	}
	else {
		throw new Error("Invalid point: point must either be array of 2 coordinates, or an L.latLng.");
	}

	return coordinate_array;
}

/**
 * Given two coordinate arrays, get their intersections.
 * 
 * @param {array<array<number>>} arr_a - Array of coordinate pairs.
 * @param {array<array<number>>} arr_b - Array of coordinate pairs.
 * @param {array<number>} ids - 2-element array whose elements are IDs for arr_a and arr_b respectively.
 *
 * @returns {Array<Array<number|Object<number, number>>>} - Array whose elements are the intersections' cooridinate-pairs if
 * ids is empty, or otherwise whose elements are arrays each of whose first element is an
 * intersection's coordinate-pair and whose second element is an object mapping each array's ID (supplied by ids) 
 * to the index of the intersection's coordinate-pair in that array.
 */
function getIntersections(arr_a, arr_b, ids = []) {
	let intersections = [];

	for (let i = 0; i < arr_a.length; i++) {
		let el_a = arr_a[i];

		for (let j = 0; j < arr_b.length; j++) {
			let el_b = arr_b[j];
			
			if (isPointCoordinates(el_a) && isPointCoordinates(el_b)) {
				if (el_a[0] === el_b[0] && el_a[1] === el_b[1]) {
					let new_intersection;

					if (ids.length === 2) {
						let identified_intersections = {};
						identified_intersections[ids[0]] = i,
						identified_intersections[ids[1]] = j,
						new_intersection = [el_a, identified_intersections];
					}
					else {
						new_intersection = el_a;
					}
				
					intersections.push(new_intersection);
				}
			}
			else {
				throw new Error("Every element of each array must be a coordinate pair array.");
			}
		}
	}

	return intersections;
}

exports.getIntersections = getIntersections;
exports.reversedCoordinates = reversedCoordinates;
exports.isPointCoordinates = isPointCoordinates;
exports.pointToCoordinateArray = pointToCoordinateArray;

},{}],49:[function(require,module,exports){
// https://d3js.org/d3-array/ v1.2.4 Copyright 2018 Mike Bostock
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(factory((global.d3 = global.d3 || {})));
}(this, (function (exports) { 'use strict';

function ascending(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

function bisector(compare) {
  if (compare.length === 1) compare = ascendingComparator(compare);
  return {
    left: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) < 0) lo = mid + 1;
        else hi = mid;
      }
      return lo;
    },
    right: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) > 0) hi = mid;
        else lo = mid + 1;
      }
      return lo;
    }
  };
}

function ascendingComparator(f) {
  return function(d, x) {
    return ascending(f(d), x);
  };
}

var ascendingBisect = bisector(ascending);
var bisectRight = ascendingBisect.right;
var bisectLeft = ascendingBisect.left;

function pairs(array, f) {
  if (f == null) f = pair;
  var i = 0, n = array.length - 1, p = array[0], pairs = new Array(n < 0 ? 0 : n);
  while (i < n) pairs[i] = f(p, p = array[++i]);
  return pairs;
}

function pair(a, b) {
  return [a, b];
}

function cross(values0, values1, reduce) {
  var n0 = values0.length,
      n1 = values1.length,
      values = new Array(n0 * n1),
      i0,
      i1,
      i,
      value0;

  if (reduce == null) reduce = pair;

  for (i0 = i = 0; i0 < n0; ++i0) {
    for (value0 = values0[i0], i1 = 0; i1 < n1; ++i1, ++i) {
      values[i] = reduce(value0, values1[i1]);
    }
  }

  return values;
}

function descending(a, b) {
  return b < a ? -1 : b > a ? 1 : b >= a ? 0 : NaN;
}

function number(x) {
  return x === null ? NaN : +x;
}

function variance(values, valueof) {
  var n = values.length,
      m = 0,
      i = -1,
      mean = 0,
      value,
      delta,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        delta = value - mean;
        mean += delta / ++m;
        sum += delta * (value - mean);
      }
    }
  }

  if (m > 1) return sum / (m - 1);
}

function deviation(array, f) {
  var v = variance(array, f);
  return v ? Math.sqrt(v) : v;
}

function extent(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  return [min, max];
}

var array = Array.prototype;

var slice = array.slice;
var map = array.map;

function constant(x) {
  return function() {
    return x;
  };
}

function identity(x) {
  return x;
}

function range(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;

  var i = -1,
      n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
      range = new Array(n);

  while (++i < n) {
    range[i] = start + i * step;
  }

  return range;
}

var e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

function ticks(start, stop, count) {
  var reverse,
      i = -1,
      n,
      ticks,
      step;

  stop = +stop, start = +start, count = +count;
  if (start === stop && count > 0) return [start];
  if (reverse = stop < start) n = start, start = stop, stop = n;
  if ((step = tickIncrement(start, stop, count)) === 0 || !isFinite(step)) return [];

  if (step > 0) {
    start = Math.ceil(start / step);
    stop = Math.floor(stop / step);
    ticks = new Array(n = Math.ceil(stop - start + 1));
    while (++i < n) ticks[i] = (start + i) * step;
  } else {
    start = Math.floor(start * step);
    stop = Math.ceil(stop * step);
    ticks = new Array(n = Math.ceil(start - stop + 1));
    while (++i < n) ticks[i] = (start - i) / step;
  }

  if (reverse) ticks.reverse();

  return ticks;
}

function tickIncrement(start, stop, count) {
  var step = (stop - start) / Math.max(0, count),
      power = Math.floor(Math.log(step) / Math.LN10),
      error = step / Math.pow(10, power);
  return power >= 0
      ? (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1) * Math.pow(10, power)
      : -Math.pow(10, -power) / (error >= e10 ? 10 : error >= e5 ? 5 : error >= e2 ? 2 : 1);
}

function tickStep(start, stop, count) {
  var step0 = Math.abs(stop - start) / Math.max(0, count),
      step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
      error = step0 / step1;
  if (error >= e10) step1 *= 10;
  else if (error >= e5) step1 *= 5;
  else if (error >= e2) step1 *= 2;
  return stop < start ? -step1 : step1;
}

function sturges(values) {
  return Math.ceil(Math.log(values.length) / Math.LN2) + 1;
}

function histogram() {
  var value = identity,
      domain = extent,
      threshold = sturges;

  function histogram(data) {
    var i,
        n = data.length,
        x,
        values = new Array(n);

    for (i = 0; i < n; ++i) {
      values[i] = value(data[i], i, data);
    }

    var xz = domain(values),
        x0 = xz[0],
        x1 = xz[1],
        tz = threshold(values, x0, x1);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      tz = tickStep(x0, x1, tz);
      tz = range(Math.ceil(x0 / tz) * tz, x1, tz); // exclusive
    }

    // Remove any thresholds outside the domain.
    var m = tz.length;
    while (tz[0] <= x0) tz.shift(), --m;
    while (tz[m - 1] > x1) tz.pop(), --m;

    var bins = new Array(m + 1),
        bin;

    // Initialize bins.
    for (i = 0; i <= m; ++i) {
      bin = bins[i] = [];
      bin.x0 = i > 0 ? tz[i - 1] : x0;
      bin.x1 = i < m ? tz[i] : x1;
    }

    // Assign data to bins by value, ignoring any outside the domain.
    for (i = 0; i < n; ++i) {
      x = values[i];
      if (x0 <= x && x <= x1) {
        bins[bisectRight(tz, x, 0, m)].push(data[i]);
      }
    }

    return bins;
  }

  histogram.value = function(_) {
    return arguments.length ? (value = typeof _ === "function" ? _ : constant(_), histogram) : value;
  };

  histogram.domain = function(_) {
    return arguments.length ? (domain = typeof _ === "function" ? _ : constant([_[0], _[1]]), histogram) : domain;
  };

  histogram.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), histogram) : threshold;
  };

  return histogram;
}

function quantile(values, p, valueof) {
  if (valueof == null) valueof = number;
  if (!(n = values.length)) return;
  if ((p = +p) <= 0 || n < 2) return +valueof(values[0], 0, values);
  if (p >= 1) return +valueof(values[n - 1], n - 1, values);
  var n,
      i = (n - 1) * p,
      i0 = Math.floor(i),
      value0 = +valueof(values[i0], i0, values),
      value1 = +valueof(values[i0 + 1], i0 + 1, values);
  return value0 + (value1 - value0) * (i - i0);
}

function freedmanDiaconis(values, min, max) {
  values = map.call(values, number).sort(ascending);
  return Math.ceil((max - min) / (2 * (quantile(values, 0.75) - quantile(values, 0.25)) * Math.pow(values.length, -1 / 3)));
}

function scott(values, min, max) {
  return Math.ceil((max - min) / (3.5 * deviation(values) * Math.pow(values.length, -1 / 3)));
}

function max(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && value > max) {
            max = value;
          }
        }
      }
    }
  }

  return max;
}

function mean(values, valueof) {
  var n = values.length,
      m = n,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) sum += value;
      else --m;
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) sum += value;
      else --m;
    }
  }

  if (m) return sum / m;
}

function median(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      numbers = [];

  if (valueof == null) {
    while (++i < n) {
      if (!isNaN(value = number(values[i]))) {
        numbers.push(value);
      }
    }
  }

  else {
    while (++i < n) {
      if (!isNaN(value = number(valueof(values[i], i, values)))) {
        numbers.push(value);
      }
    }
  }

  return quantile(numbers.sort(ascending), 0.5);
}

function merge(arrays) {
  var n = arrays.length,
      m,
      i = -1,
      j = 0,
      merged,
      array;

  while (++i < n) j += arrays[i].length;
  merged = new Array(j);

  while (--n >= 0) {
    array = arrays[n];
    m = array.length;
    while (--m >= 0) {
      merged[--j] = array[m];
    }
  }

  return merged;
}

function min(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null && min > value) {
            min = value;
          }
        }
      }
    }
  }

  return min;
}

function permute(array, indexes) {
  var i = indexes.length, permutes = new Array(i);
  while (i--) permutes[i] = array[indexes[i]];
  return permutes;
}

function scan(values, compare) {
  if (!(n = values.length)) return;
  var n,
      i = 0,
      j = 0,
      xi,
      xj = values[j];

  if (compare == null) compare = ascending;

  while (++i < n) {
    if (compare(xi = values[i], xj) < 0 || compare(xj, xj) !== 0) {
      xj = xi, j = i;
    }
  }

  if (compare(xj, xj) === 0) return j;
}

function shuffle(array, i0, i1) {
  var m = (i1 == null ? array.length : i1) - (i0 = i0 == null ? 0 : +i0),
      t,
      i;

  while (m) {
    i = Math.random() * m-- | 0;
    t = array[m + i0];
    array[m + i0] = array[i + i0];
    array[i + i0] = t;
  }

  return array;
}

function sum(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      sum = 0;

  if (valueof == null) {
    while (++i < n) {
      if (value = +values[i]) sum += value; // Note: zero and null are equivalent.
    }
  }

  else {
    while (++i < n) {
      if (value = +valueof(values[i], i, values)) sum += value;
    }
  }

  return sum;
}

function transpose(matrix) {
  if (!(n = matrix.length)) return [];
  for (var i = -1, m = min(matrix, length), transpose = new Array(m); ++i < m;) {
    for (var j = -1, n, row = transpose[i] = new Array(n); ++j < n;) {
      row[j] = matrix[j][i];
    }
  }
  return transpose;
}

function length(d) {
  return d.length;
}

function zip() {
  return transpose(arguments);
}

exports.bisect = bisectRight;
exports.bisectRight = bisectRight;
exports.bisectLeft = bisectLeft;
exports.ascending = ascending;
exports.bisector = bisector;
exports.cross = cross;
exports.descending = descending;
exports.deviation = deviation;
exports.extent = extent;
exports.histogram = histogram;
exports.thresholdFreedmanDiaconis = freedmanDiaconis;
exports.thresholdScott = scott;
exports.thresholdSturges = sturges;
exports.max = max;
exports.mean = mean;
exports.median = median;
exports.merge = merge;
exports.min = min;
exports.pairs = pairs;
exports.permute = permute;
exports.quantile = quantile;
exports.range = range;
exports.scan = scan;
exports.shuffle = shuffle;
exports.sum = sum;
exports.ticks = ticks;
exports.tickIncrement = tickIncrement;
exports.tickStep = tickStep;
exports.transpose = transpose;
exports.variance = variance;
exports.zip = zip;

Object.defineProperty(exports, '__esModule', { value: true });

})));

},{}],50:[function(require,module,exports){
// https://d3js.org/d3-geo/ Version 1.7.1. Copyright 2017 Mike Bostock.
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('d3-array')) :
	typeof define === 'function' && define.amd ? define(['exports', 'd3-array'], factory) :
	(factory((global.d3 = global.d3 || {}),global.d3));
}(this, (function (exports,d3Array) { 'use strict';

// Adds floating point numbers with twice the normal precision.
// Reference: J. R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and
// Fast Robust Geometric Predicates, Discrete & Computational Geometry 18(3)
// 305–363 (1997).
// Code adapted from GeographicLib by Charles F. F. Karney,
// http://geographiclib.sourceforge.net/

var adder = function() {
  return new Adder;
};

function Adder() {
  this.reset();
}

Adder.prototype = {
  constructor: Adder,
  reset: function() {
    this.s = // rounded value
    this.t = 0; // exact error
  },
  add: function(y) {
    add(temp, y, this.t);
    add(this, temp.s, this.s);
    if (this.s) this.t += temp.t;
    else this.s = temp.t;
  },
  valueOf: function() {
    return this.s;
  }
};

var temp = new Adder;

function add(adder, a, b) {
  var x = adder.s = a + b,
      bv = x - a,
      av = x - bv;
  adder.t = (a - av) + (b - bv);
}

var epsilon = 1e-6;
var epsilon2 = 1e-12;
var pi = Math.PI;
var halfPi = pi / 2;
var quarterPi = pi / 4;
var tau = pi * 2;

var degrees = 180 / pi;
var radians = pi / 180;

var abs = Math.abs;
var atan = Math.atan;
var atan2 = Math.atan2;
var cos = Math.cos;
var ceil = Math.ceil;
var exp = Math.exp;

var log = Math.log;
var pow = Math.pow;
var sin = Math.sin;
var sign = Math.sign || function(x) { return x > 0 ? 1 : x < 0 ? -1 : 0; };
var sqrt = Math.sqrt;
var tan = Math.tan;

function acos(x) {
  return x > 1 ? 0 : x < -1 ? pi : Math.acos(x);
}

function asin(x) {
  return x > 1 ? halfPi : x < -1 ? -halfPi : Math.asin(x);
}

function haversin(x) {
  return (x = sin(x / 2)) * x;
}

function noop() {}

function streamGeometry(geometry, stream) {
  if (geometry && streamGeometryType.hasOwnProperty(geometry.type)) {
    streamGeometryType[geometry.type](geometry, stream);
  }
}

var streamObjectType = {
  Feature: function(object, stream) {
    streamGeometry(object.geometry, stream);
  },
  FeatureCollection: function(object, stream) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) streamGeometry(features[i].geometry, stream);
  }
};

var streamGeometryType = {
  Sphere: function(object, stream) {
    stream.sphere();
  },
  Point: function(object, stream) {
    object = object.coordinates;
    stream.point(object[0], object[1], object[2]);
  },
  MultiPoint: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) object = coordinates[i], stream.point(object[0], object[1], object[2]);
  },
  LineString: function(object, stream) {
    streamLine(object.coordinates, stream, 0);
  },
  MultiLineString: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamLine(coordinates[i], stream, 0);
  },
  Polygon: function(object, stream) {
    streamPolygon(object.coordinates, stream);
  },
  MultiPolygon: function(object, stream) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) streamPolygon(coordinates[i], stream);
  },
  GeometryCollection: function(object, stream) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) streamGeometry(geometries[i], stream);
  }
};

function streamLine(coordinates, stream, closed) {
  var i = -1, n = coordinates.length - closed, coordinate;
  stream.lineStart();
  while (++i < n) coordinate = coordinates[i], stream.point(coordinate[0], coordinate[1], coordinate[2]);
  stream.lineEnd();
}

function streamPolygon(coordinates, stream) {
  var i = -1, n = coordinates.length;
  stream.polygonStart();
  while (++i < n) streamLine(coordinates[i], stream, 1);
  stream.polygonEnd();
}

var geoStream = function(object, stream) {
  if (object && streamObjectType.hasOwnProperty(object.type)) {
    streamObjectType[object.type](object, stream);
  } else {
    streamGeometry(object, stream);
  }
};

var areaRingSum = adder();

var areaSum = adder();
var lambda00;
var phi00;
var lambda0;
var cosPhi0;
var sinPhi0;

var areaStream = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaRingSum.reset();
    areaStream.lineStart = areaRingStart;
    areaStream.lineEnd = areaRingEnd;
  },
  polygonEnd: function() {
    var areaRing = +areaRingSum;
    areaSum.add(areaRing < 0 ? tau + areaRing : areaRing);
    this.lineStart = this.lineEnd = this.point = noop;
  },
  sphere: function() {
    areaSum.add(tau);
  }
};

function areaRingStart() {
  areaStream.point = areaPointFirst;
}

function areaRingEnd() {
  areaPoint(lambda00, phi00);
}

function areaPointFirst(lambda, phi) {
  areaStream.point = areaPoint;
  lambda00 = lambda, phi00 = phi;
  lambda *= radians, phi *= radians;
  lambda0 = lambda, cosPhi0 = cos(phi = phi / 2 + quarterPi), sinPhi0 = sin(phi);
}

function areaPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  phi = phi / 2 + quarterPi; // half the angular distance from south pole

  // Spherical excess E for a spherical triangle with vertices: south pole,
  // previous point, current point.  Uses a formula derived from Cagnoli’s
  // theorem.  See Todhunter, Spherical Trig. (1871), Sec. 103, Eq. (2).
  var dLambda = lambda - lambda0,
      sdLambda = dLambda >= 0 ? 1 : -1,
      adLambda = sdLambda * dLambda,
      cosPhi = cos(phi),
      sinPhi = sin(phi),
      k = sinPhi0 * sinPhi,
      u = cosPhi0 * cosPhi + k * cos(adLambda),
      v = k * sdLambda * sin(adLambda);
  areaRingSum.add(atan2(v, u));

  // Advance the previous points.
  lambda0 = lambda, cosPhi0 = cosPhi, sinPhi0 = sinPhi;
}

var area = function(object) {
  areaSum.reset();
  geoStream(object, areaStream);
  return areaSum * 2;
};

function spherical(cartesian) {
  return [atan2(cartesian[1], cartesian[0]), asin(cartesian[2])];
}

function cartesian(spherical) {
  var lambda = spherical[0], phi = spherical[1], cosPhi = cos(phi);
  return [cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi)];
}

function cartesianDot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function cartesianCross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

// TODO return a
function cartesianAddInPlace(a, b) {
  a[0] += b[0], a[1] += b[1], a[2] += b[2];
}

function cartesianScale(vector, k) {
  return [vector[0] * k, vector[1] * k, vector[2] * k];
}

// TODO return d
function cartesianNormalizeInPlace(d) {
  var l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  d[0] /= l, d[1] /= l, d[2] /= l;
}

var lambda0$1;
var phi0;
var lambda1;
var phi1;
var lambda2;
var lambda00$1;
var phi00$1;
var p0;
var deltaSum = adder();
var ranges;
var range$1;

var boundsStream = {
  point: boundsPoint,
  lineStart: boundsLineStart,
  lineEnd: boundsLineEnd,
  polygonStart: function() {
    boundsStream.point = boundsRingPoint;
    boundsStream.lineStart = boundsRingStart;
    boundsStream.lineEnd = boundsRingEnd;
    deltaSum.reset();
    areaStream.polygonStart();
  },
  polygonEnd: function() {
    areaStream.polygonEnd();
    boundsStream.point = boundsPoint;
    boundsStream.lineStart = boundsLineStart;
    boundsStream.lineEnd = boundsLineEnd;
    if (areaRingSum < 0) lambda0$1 = -(lambda1 = 180), phi0 = -(phi1 = 90);
    else if (deltaSum > epsilon) phi1 = 90;
    else if (deltaSum < -epsilon) phi0 = -90;
    range$1[0] = lambda0$1, range$1[1] = lambda1;
  }
};

function boundsPoint(lambda, phi) {
  ranges.push(range$1 = [lambda0$1 = lambda, lambda1 = lambda]);
  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
}

function linePoint(lambda, phi) {
  var p = cartesian([lambda * radians, phi * radians]);
  if (p0) {
    var normal = cartesianCross(p0, p),
        equatorial = [normal[1], -normal[0], 0],
        inflection = cartesianCross(equatorial, normal);
    cartesianNormalizeInPlace(inflection);
    inflection = spherical(inflection);
    var delta = lambda - lambda2,
        sign$$1 = delta > 0 ? 1 : -1,
        lambdai = inflection[0] * degrees * sign$$1,
        phii,
        antimeridian = abs(delta) > 180;
    if (antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
      phii = inflection[1] * degrees;
      if (phii > phi1) phi1 = phii;
    } else if (lambdai = (lambdai + 360) % 360 - 180, antimeridian ^ (sign$$1 * lambda2 < lambdai && lambdai < sign$$1 * lambda)) {
      phii = -inflection[1] * degrees;
      if (phii < phi0) phi0 = phii;
    } else {
      if (phi < phi0) phi0 = phi;
      if (phi > phi1) phi1 = phi;
    }
    if (antimeridian) {
      if (lambda < lambda2) {
        if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
      } else {
        if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
      }
    } else {
      if (lambda1 >= lambda0$1) {
        if (lambda < lambda0$1) lambda0$1 = lambda;
        if (lambda > lambda1) lambda1 = lambda;
      } else {
        if (lambda > lambda2) {
          if (angle(lambda0$1, lambda) > angle(lambda0$1, lambda1)) lambda1 = lambda;
        } else {
          if (angle(lambda, lambda1) > angle(lambda0$1, lambda1)) lambda0$1 = lambda;
        }
      }
    }
  } else {
    ranges.push(range$1 = [lambda0$1 = lambda, lambda1 = lambda]);
  }
  if (phi < phi0) phi0 = phi;
  if (phi > phi1) phi1 = phi;
  p0 = p, lambda2 = lambda;
}

function boundsLineStart() {
  boundsStream.point = linePoint;
}

function boundsLineEnd() {
  range$1[0] = lambda0$1, range$1[1] = lambda1;
  boundsStream.point = boundsPoint;
  p0 = null;
}

function boundsRingPoint(lambda, phi) {
  if (p0) {
    var delta = lambda - lambda2;
    deltaSum.add(abs(delta) > 180 ? delta + (delta > 0 ? 360 : -360) : delta);
  } else {
    lambda00$1 = lambda, phi00$1 = phi;
  }
  areaStream.point(lambda, phi);
  linePoint(lambda, phi);
}

function boundsRingStart() {
  areaStream.lineStart();
}

function boundsRingEnd() {
  boundsRingPoint(lambda00$1, phi00$1);
  areaStream.lineEnd();
  if (abs(deltaSum) > epsilon) lambda0$1 = -(lambda1 = 180);
  range$1[0] = lambda0$1, range$1[1] = lambda1;
  p0 = null;
}

// Finds the left-right distance between two longitudes.
// This is almost the same as (lambda1 - lambda0 + 360°) % 360°, except that we want
// the distance between ±180° to be 360°.
function angle(lambda0, lambda1) {
  return (lambda1 -= lambda0) < 0 ? lambda1 + 360 : lambda1;
}

function rangeCompare(a, b) {
  return a[0] - b[0];
}

function rangeContains(range$$1, x) {
  return range$$1[0] <= range$$1[1] ? range$$1[0] <= x && x <= range$$1[1] : x < range$$1[0] || range$$1[1] < x;
}

var bounds = function(feature) {
  var i, n, a, b, merged, deltaMax, delta;

  phi1 = lambda1 = -(lambda0$1 = phi0 = Infinity);
  ranges = [];
  geoStream(feature, boundsStream);

  // First, sort ranges by their minimum longitudes.
  if (n = ranges.length) {
    ranges.sort(rangeCompare);

    // Then, merge any ranges that overlap.
    for (i = 1, a = ranges[0], merged = [a]; i < n; ++i) {
      b = ranges[i];
      if (rangeContains(a, b[0]) || rangeContains(a, b[1])) {
        if (angle(a[0], b[1]) > angle(a[0], a[1])) a[1] = b[1];
        if (angle(b[0], a[1]) > angle(a[0], a[1])) a[0] = b[0];
      } else {
        merged.push(a = b);
      }
    }

    // Finally, find the largest gap between the merged ranges.
    // The final bounding box will be the inverse of this gap.
    for (deltaMax = -Infinity, n = merged.length - 1, i = 0, a = merged[n]; i <= n; a = b, ++i) {
      b = merged[i];
      if ((delta = angle(a[1], b[0])) > deltaMax) deltaMax = delta, lambda0$1 = b[0], lambda1 = a[1];
    }
  }

  ranges = range$1 = null;

  return lambda0$1 === Infinity || phi0 === Infinity
      ? [[NaN, NaN], [NaN, NaN]]
      : [[lambda0$1, phi0], [lambda1, phi1]];
};

var W0;
var W1;
var X0;
var Y0;
var Z0;
var X1;
var Y1;
var Z1;
var X2;
var Y2;
var Z2;
var lambda00$2;
var phi00$2;
var x0;
var y0;
var z0; // previous point

var centroidStream = {
  sphere: noop,
  point: centroidPoint,
  lineStart: centroidLineStart,
  lineEnd: centroidLineEnd,
  polygonStart: function() {
    centroidStream.lineStart = centroidRingStart;
    centroidStream.lineEnd = centroidRingEnd;
  },
  polygonEnd: function() {
    centroidStream.lineStart = centroidLineStart;
    centroidStream.lineEnd = centroidLineEnd;
  }
};

// Arithmetic mean of Cartesian vectors.
function centroidPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  centroidPointCartesian(cosPhi * cos(lambda), cosPhi * sin(lambda), sin(phi));
}

function centroidPointCartesian(x, y, z) {
  ++W0;
  X0 += (x - X0) / W0;
  Y0 += (y - Y0) / W0;
  Z0 += (z - Z0) / W0;
}

function centroidLineStart() {
  centroidStream.point = centroidLinePointFirst;
}

function centroidLinePointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidStream.point = centroidLinePoint;
  centroidPointCartesian(x0, y0, z0);
}

function centroidLinePoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      w = atan2(sqrt((w = y0 * z - z0 * y) * w + (w = z0 * x - x0 * z) * w + (w = x0 * y - y0 * x) * w), x0 * x + y0 * y + z0 * z);
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

function centroidLineEnd() {
  centroidStream.point = centroidPoint;
}

// See J. E. Brock, The Inertia Tensor for a Spherical Triangle,
// J. Applied Mechanics 42, 239 (1975).
function centroidRingStart() {
  centroidStream.point = centroidRingPointFirst;
}

function centroidRingEnd() {
  centroidRingPoint(lambda00$2, phi00$2);
  centroidStream.point = centroidPoint;
}

function centroidRingPointFirst(lambda, phi) {
  lambda00$2 = lambda, phi00$2 = phi;
  lambda *= radians, phi *= radians;
  centroidStream.point = centroidRingPoint;
  var cosPhi = cos(phi);
  x0 = cosPhi * cos(lambda);
  y0 = cosPhi * sin(lambda);
  z0 = sin(phi);
  centroidPointCartesian(x0, y0, z0);
}

function centroidRingPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var cosPhi = cos(phi),
      x = cosPhi * cos(lambda),
      y = cosPhi * sin(lambda),
      z = sin(phi),
      cx = y0 * z - z0 * y,
      cy = z0 * x - x0 * z,
      cz = x0 * y - y0 * x,
      m = sqrt(cx * cx + cy * cy + cz * cz),
      w = asin(m), // line weight = angle
      v = m && -w / m; // area weight multiplier
  X2 += v * cx;
  Y2 += v * cy;
  Z2 += v * cz;
  W1 += w;
  X1 += w * (x0 + (x0 = x));
  Y1 += w * (y0 + (y0 = y));
  Z1 += w * (z0 + (z0 = z));
  centroidPointCartesian(x0, y0, z0);
}

var centroid = function(object) {
  W0 = W1 =
  X0 = Y0 = Z0 =
  X1 = Y1 = Z1 =
  X2 = Y2 = Z2 = 0;
  geoStream(object, centroidStream);

  var x = X2,
      y = Y2,
      z = Z2,
      m = x * x + y * y + z * z;

  // If the area-weighted ccentroid is undefined, fall back to length-weighted ccentroid.
  if (m < epsilon2) {
    x = X1, y = Y1, z = Z1;
    // If the feature has zero length, fall back to arithmetic mean of point vectors.
    if (W1 < epsilon) x = X0, y = Y0, z = Z0;
    m = x * x + y * y + z * z;
    // If the feature still has an undefined ccentroid, then return.
    if (m < epsilon2) return [NaN, NaN];
  }

  return [atan2(y, x) * degrees, asin(z / sqrt(m)) * degrees];
};

var constant = function(x) {
  return function() {
    return x;
  };
};

var compose = function(a, b) {

  function compose(x, y) {
    return x = a(x, y), b(x[0], x[1]);
  }

  if (a.invert && b.invert) compose.invert = function(x, y) {
    return x = b.invert(x, y), x && a.invert(x[0], x[1]);
  };

  return compose;
};

function rotationIdentity(lambda, phi) {
  return [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
}

rotationIdentity.invert = rotationIdentity;

function rotateRadians(deltaLambda, deltaPhi, deltaGamma) {
  return (deltaLambda %= tau) ? (deltaPhi || deltaGamma ? compose(rotationLambda(deltaLambda), rotationPhiGamma(deltaPhi, deltaGamma))
    : rotationLambda(deltaLambda))
    : (deltaPhi || deltaGamma ? rotationPhiGamma(deltaPhi, deltaGamma)
    : rotationIdentity);
}

function forwardRotationLambda(deltaLambda) {
  return function(lambda, phi) {
    return lambda += deltaLambda, [lambda > pi ? lambda - tau : lambda < -pi ? lambda + tau : lambda, phi];
  };
}

function rotationLambda(deltaLambda) {
  var rotation = forwardRotationLambda(deltaLambda);
  rotation.invert = forwardRotationLambda(-deltaLambda);
  return rotation;
}

function rotationPhiGamma(deltaPhi, deltaGamma) {
  var cosDeltaPhi = cos(deltaPhi),
      sinDeltaPhi = sin(deltaPhi),
      cosDeltaGamma = cos(deltaGamma),
      sinDeltaGamma = sin(deltaGamma);

  function rotation(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaPhi + x * sinDeltaPhi;
    return [
      atan2(y * cosDeltaGamma - k * sinDeltaGamma, x * cosDeltaPhi - z * sinDeltaPhi),
      asin(k * cosDeltaGamma + y * sinDeltaGamma)
    ];
  }

  rotation.invert = function(lambda, phi) {
    var cosPhi = cos(phi),
        x = cos(lambda) * cosPhi,
        y = sin(lambda) * cosPhi,
        z = sin(phi),
        k = z * cosDeltaGamma - y * sinDeltaGamma;
    return [
      atan2(y * cosDeltaGamma + z * sinDeltaGamma, x * cosDeltaPhi + k * sinDeltaPhi),
      asin(k * cosDeltaPhi - x * sinDeltaPhi)
    ];
  };

  return rotation;
}

var rotation = function(rotate) {
  rotate = rotateRadians(rotate[0] * radians, rotate[1] * radians, rotate.length > 2 ? rotate[2] * radians : 0);

  function forward(coordinates) {
    coordinates = rotate(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  }

  forward.invert = function(coordinates) {
    coordinates = rotate.invert(coordinates[0] * radians, coordinates[1] * radians);
    return coordinates[0] *= degrees, coordinates[1] *= degrees, coordinates;
  };

  return forward;
};

// Generates a circle centered at [0°, 0°], with a given radius and precision.
function circleStream(stream, radius, delta, direction, t0, t1) {
  if (!delta) return;
  var cosRadius = cos(radius),
      sinRadius = sin(radius),
      step = direction * delta;
  if (t0 == null) {
    t0 = radius + direction * tau;
    t1 = radius - step / 2;
  } else {
    t0 = circleRadius(cosRadius, t0);
    t1 = circleRadius(cosRadius, t1);
    if (direction > 0 ? t0 < t1 : t0 > t1) t0 += direction * tau;
  }
  for (var point, t = t0; direction > 0 ? t > t1 : t < t1; t -= step) {
    point = spherical([cosRadius, -sinRadius * cos(t), -sinRadius * sin(t)]);
    stream.point(point[0], point[1]);
  }
}

// Returns the signed angle of a cartesian point relative to [cosRadius, 0, 0].
function circleRadius(cosRadius, point) {
  point = cartesian(point), point[0] -= cosRadius;
  cartesianNormalizeInPlace(point);
  var radius = acos(-point[1]);
  return ((-point[2] < 0 ? -radius : radius) + tau - epsilon) % tau;
}

var circle = function() {
  var center = constant([0, 0]),
      radius = constant(90),
      precision = constant(6),
      ring,
      rotate,
      stream = {point: point};

  function point(x, y) {
    ring.push(x = rotate(x, y));
    x[0] *= degrees, x[1] *= degrees;
  }

  function circle() {
    var c = center.apply(this, arguments),
        r = radius.apply(this, arguments) * radians,
        p = precision.apply(this, arguments) * radians;
    ring = [];
    rotate = rotateRadians(-c[0] * radians, -c[1] * radians, 0).invert;
    circleStream(stream, r, p, 1);
    c = {type: "Polygon", coordinates: [ring]};
    ring = rotate = null;
    return c;
  }

  circle.center = function(_) {
    return arguments.length ? (center = typeof _ === "function" ? _ : constant([+_[0], +_[1]]), circle) : center;
  };

  circle.radius = function(_) {
    return arguments.length ? (radius = typeof _ === "function" ? _ : constant(+_), circle) : radius;
  };

  circle.precision = function(_) {
    return arguments.length ? (precision = typeof _ === "function" ? _ : constant(+_), circle) : precision;
  };

  return circle;
};

var clipBuffer = function() {
  var lines = [],
      line;
  return {
    point: function(x, y) {
      line.push([x, y]);
    },
    lineStart: function() {
      lines.push(line = []);
    },
    lineEnd: noop,
    rejoin: function() {
      if (lines.length > 1) lines.push(lines.pop().concat(lines.shift()));
    },
    result: function() {
      var result = lines;
      lines = [];
      line = null;
      return result;
    }
  };
};

var clipLine = function(a, b, x0, y0, x1, y1) {
  var ax = a[0],
      ay = a[1],
      bx = b[0],
      by = b[1],
      t0 = 0,
      t1 = 1,
      dx = bx - ax,
      dy = by - ay,
      r;

  r = x0 - ax;
  if (!dx && r > 0) return;
  r /= dx;
  if (dx < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dx > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = x1 - ax;
  if (!dx && r < 0) return;
  r /= dx;
  if (dx < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dx > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  r = y0 - ay;
  if (!dy && r > 0) return;
  r /= dy;
  if (dy < 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  } else if (dy > 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  }

  r = y1 - ay;
  if (!dy && r < 0) return;
  r /= dy;
  if (dy < 0) {
    if (r > t1) return;
    if (r > t0) t0 = r;
  } else if (dy > 0) {
    if (r < t0) return;
    if (r < t1) t1 = r;
  }

  if (t0 > 0) a[0] = ax + t0 * dx, a[1] = ay + t0 * dy;
  if (t1 < 1) b[0] = ax + t1 * dx, b[1] = ay + t1 * dy;
  return true;
};

var pointEqual = function(a, b) {
  return abs(a[0] - b[0]) < epsilon && abs(a[1] - b[1]) < epsilon;
};

function Intersection(point, points, other, entry) {
  this.x = point;
  this.z = points;
  this.o = other; // another intersection
  this.e = entry; // is an entry?
  this.v = false; // visited
  this.n = this.p = null; // next & previous
}

// A generalized polygon clipping algorithm: given a polygon that has been cut
// into its visible line segments, and rejoins the segments by interpolating
// along the clip edge.
var clipPolygon = function(segments, compareIntersection, startInside, interpolate, stream) {
  var subject = [],
      clip = [],
      i,
      n;

  segments.forEach(function(segment) {
    if ((n = segment.length - 1) <= 0) return;
    var n, p0 = segment[0], p1 = segment[n], x;

    // If the first and last points of a segment are coincident, then treat as a
    // closed ring. TODO if all rings are closed, then the winding order of the
    // exterior ring should be checked.
    if (pointEqual(p0, p1)) {
      stream.lineStart();
      for (i = 0; i < n; ++i) stream.point((p0 = segment[i])[0], p0[1]);
      stream.lineEnd();
      return;
    }

    subject.push(x = new Intersection(p0, segment, null, true));
    clip.push(x.o = new Intersection(p0, null, x, false));
    subject.push(x = new Intersection(p1, segment, null, false));
    clip.push(x.o = new Intersection(p1, null, x, true));
  });

  if (!subject.length) return;

  clip.sort(compareIntersection);
  link(subject);
  link(clip);

  for (i = 0, n = clip.length; i < n; ++i) {
    clip[i].e = startInside = !startInside;
  }

  var start = subject[0],
      points,
      point;

  while (1) {
    // Find first unvisited intersection.
    var current = start,
        isSubject = true;
    while (current.v) if ((current = current.n) === start) return;
    points = current.z;
    stream.lineStart();
    do {
      current.v = current.o.v = true;
      if (current.e) {
        if (isSubject) {
          for (i = 0, n = points.length; i < n; ++i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.n.x, 1, stream);
        }
        current = current.n;
      } else {
        if (isSubject) {
          points = current.p.z;
          for (i = points.length - 1; i >= 0; --i) stream.point((point = points[i])[0], point[1]);
        } else {
          interpolate(current.x, current.p.x, -1, stream);
        }
        current = current.p;
      }
      current = current.o;
      points = current.z;
      isSubject = !isSubject;
    } while (!current.v);
    stream.lineEnd();
  }
};

function link(array) {
  if (!(n = array.length)) return;
  var n,
      i = 0,
      a = array[0],
      b;
  while (++i < n) {
    a.n = b = array[i];
    b.p = a;
    a = b;
  }
  a.n = b = array[0];
  b.p = a;
}

var clipMax = 1e9;
var clipMin = -clipMax;

// TODO Use d3-polygon’s polygonContains here for the ring check?
// TODO Eliminate duplicate buffering in clipBuffer and polygon.push?

function clipExtent(x0, y0, x1, y1) {

  function visible(x, y) {
    return x0 <= x && x <= x1 && y0 <= y && y <= y1;
  }

  function interpolate(from, to, direction, stream) {
    var a = 0, a1 = 0;
    if (from == null
        || (a = corner(from, direction)) !== (a1 = corner(to, direction))
        || comparePoint(from, to) < 0 ^ direction > 0) {
      do stream.point(a === 0 || a === 3 ? x0 : x1, a > 1 ? y1 : y0);
      while ((a = (a + direction + 4) % 4) !== a1);
    } else {
      stream.point(to[0], to[1]);
    }
  }

  function corner(p, direction) {
    return abs(p[0] - x0) < epsilon ? direction > 0 ? 0 : 3
        : abs(p[0] - x1) < epsilon ? direction > 0 ? 2 : 1
        : abs(p[1] - y0) < epsilon ? direction > 0 ? 1 : 0
        : direction > 0 ? 3 : 2; // abs(p[1] - y1) < epsilon
  }

  function compareIntersection(a, b) {
    return comparePoint(a.x, b.x);
  }

  function comparePoint(a, b) {
    var ca = corner(a, 1),
        cb = corner(b, 1);
    return ca !== cb ? ca - cb
        : ca === 0 ? b[1] - a[1]
        : ca === 1 ? a[0] - b[0]
        : ca === 2 ? a[1] - b[1]
        : b[0] - a[0];
  }

  return function(stream) {
    var activeStream = stream,
        bufferStream = clipBuffer(),
        segments,
        polygon,
        ring,
        x__, y__, v__, // first point
        x_, y_, v_, // previous point
        first,
        clean;

    var clipStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: polygonStart,
      polygonEnd: polygonEnd
    };

    function point(x, y) {
      if (visible(x, y)) activeStream.point(x, y);
    }

    function polygonInside() {
      var winding = 0;

      for (var i = 0, n = polygon.length; i < n; ++i) {
        for (var ring = polygon[i], j = 1, m = ring.length, point = ring[0], a0, a1, b0 = point[0], b1 = point[1]; j < m; ++j) {
          a0 = b0, a1 = b1, point = ring[j], b0 = point[0], b1 = point[1];
          if (a1 <= y1) { if (b1 > y1 && (b0 - a0) * (y1 - a1) > (b1 - a1) * (x0 - a0)) ++winding; }
          else { if (b1 <= y1 && (b0 - a0) * (y1 - a1) < (b1 - a1) * (x0 - a0)) --winding; }
        }
      }

      return winding;
    }

    // Buffer geometry within a polygon and then clip it en masse.
    function polygonStart() {
      activeStream = bufferStream, segments = [], polygon = [], clean = true;
    }

    function polygonEnd() {
      var startInside = polygonInside(),
          cleanInside = clean && startInside,
          visible = (segments = d3Array.merge(segments)).length;
      if (cleanInside || visible) {
        stream.polygonStart();
        if (cleanInside) {
          stream.lineStart();
          interpolate(null, null, 1, stream);
          stream.lineEnd();
        }
        if (visible) {
          clipPolygon(segments, compareIntersection, startInside, interpolate, stream);
        }
        stream.polygonEnd();
      }
      activeStream = stream, segments = polygon = ring = null;
    }

    function lineStart() {
      clipStream.point = linePoint;
      if (polygon) polygon.push(ring = []);
      first = true;
      v_ = false;
      x_ = y_ = NaN;
    }

    // TODO rather than special-case polygons, simply handle them separately.
    // Ideally, coincident intersection points should be jittered to avoid
    // clipping issues.
    function lineEnd() {
      if (segments) {
        linePoint(x__, y__);
        if (v__ && v_) bufferStream.rejoin();
        segments.push(bufferStream.result());
      }
      clipStream.point = point;
      if (v_) activeStream.lineEnd();
    }

    function linePoint(x, y) {
      var v = visible(x, y);
      if (polygon) ring.push([x, y]);
      if (first) {
        x__ = x, y__ = y, v__ = v;
        first = false;
        if (v) {
          activeStream.lineStart();
          activeStream.point(x, y);
        }
      } else {
        if (v && v_) activeStream.point(x, y);
        else {
          var a = [x_ = Math.max(clipMin, Math.min(clipMax, x_)), y_ = Math.max(clipMin, Math.min(clipMax, y_))],
              b = [x = Math.max(clipMin, Math.min(clipMax, x)), y = Math.max(clipMin, Math.min(clipMax, y))];
          if (clipLine(a, b, x0, y0, x1, y1)) {
            if (!v_) {
              activeStream.lineStart();
              activeStream.point(a[0], a[1]);
            }
            activeStream.point(b[0], b[1]);
            if (!v) activeStream.lineEnd();
            clean = false;
          } else if (v) {
            activeStream.lineStart();
            activeStream.point(x, y);
            clean = false;
          }
        }
      }
      x_ = x, y_ = y, v_ = v;
    }

    return clipStream;
  };
}

var extent = function() {
  var x0 = 0,
      y0 = 0,
      x1 = 960,
      y1 = 500,
      cache,
      cacheStream,
      clip;

  return clip = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = clipExtent(x0, y0, x1, y1)(cacheStream = stream);
    },
    extent: function(_) {
      return arguments.length ? (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1], cache = cacheStream = null, clip) : [[x0, y0], [x1, y1]];
    }
  };
};

var sum = adder();

var polygonContains = function(polygon, point) {
  var lambda = point[0],
      phi = point[1],
      normal = [sin(lambda), -cos(lambda), 0],
      angle = 0,
      winding = 0;

  sum.reset();

  for (var i = 0, n = polygon.length; i < n; ++i) {
    if (!(m = (ring = polygon[i]).length)) continue;
    var ring,
        m,
        point0 = ring[m - 1],
        lambda0 = point0[0],
        phi0 = point0[1] / 2 + quarterPi,
        sinPhi0 = sin(phi0),
        cosPhi0 = cos(phi0);

    for (var j = 0; j < m; ++j, lambda0 = lambda1, sinPhi0 = sinPhi1, cosPhi0 = cosPhi1, point0 = point1) {
      var point1 = ring[j],
          lambda1 = point1[0],
          phi1 = point1[1] / 2 + quarterPi,
          sinPhi1 = sin(phi1),
          cosPhi1 = cos(phi1),
          delta = lambda1 - lambda0,
          sign$$1 = delta >= 0 ? 1 : -1,
          absDelta = sign$$1 * delta,
          antimeridian = absDelta > pi,
          k = sinPhi0 * sinPhi1;

      sum.add(atan2(k * sign$$1 * sin(absDelta), cosPhi0 * cosPhi1 + k * cos(absDelta)));
      angle += antimeridian ? delta + sign$$1 * tau : delta;

      // Are the longitudes either side of the point’s meridian (lambda),
      // and are the latitudes smaller than the parallel (phi)?
      if (antimeridian ^ lambda0 >= lambda ^ lambda1 >= lambda) {
        var arc = cartesianCross(cartesian(point0), cartesian(point1));
        cartesianNormalizeInPlace(arc);
        var intersection = cartesianCross(normal, arc);
        cartesianNormalizeInPlace(intersection);
        var phiArc = (antimeridian ^ delta >= 0 ? -1 : 1) * asin(intersection[2]);
        if (phi > phiArc || phi === phiArc && (arc[0] || arc[1])) {
          winding += antimeridian ^ delta >= 0 ? 1 : -1;
        }
      }
    }
  }

  // First, determine whether the South pole is inside or outside:
  //
  // It is inside if:
  // * the polygon winds around it in a clockwise direction.
  // * the polygon does not (cumulatively) wind around it, but has a negative
  //   (counter-clockwise) area.
  //
  // Second, count the (signed) number of times a segment crosses a lambda
  // from the point to the South pole.  If it is zero, then the point is the
  // same side as the South pole.

  return (angle < -epsilon || angle < epsilon && sum < -epsilon) ^ (winding & 1);
};

var lengthSum = adder();
var lambda0$2;
var sinPhi0$1;
var cosPhi0$1;

var lengthStream = {
  sphere: noop,
  point: noop,
  lineStart: lengthLineStart,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop
};

function lengthLineStart() {
  lengthStream.point = lengthPointFirst;
  lengthStream.lineEnd = lengthLineEnd;
}

function lengthLineEnd() {
  lengthStream.point = lengthStream.lineEnd = noop;
}

function lengthPointFirst(lambda, phi) {
  lambda *= radians, phi *= radians;
  lambda0$2 = lambda, sinPhi0$1 = sin(phi), cosPhi0$1 = cos(phi);
  lengthStream.point = lengthPoint;
}

function lengthPoint(lambda, phi) {
  lambda *= radians, phi *= radians;
  var sinPhi = sin(phi),
      cosPhi = cos(phi),
      delta = abs(lambda - lambda0$2),
      cosDelta = cos(delta),
      sinDelta = sin(delta),
      x = cosPhi * sinDelta,
      y = cosPhi0$1 * sinPhi - sinPhi0$1 * cosPhi * cosDelta,
      z = sinPhi0$1 * sinPhi + cosPhi0$1 * cosPhi * cosDelta;
  lengthSum.add(atan2(sqrt(x * x + y * y), z));
  lambda0$2 = lambda, sinPhi0$1 = sinPhi, cosPhi0$1 = cosPhi;
}

var length = function(object) {
  lengthSum.reset();
  geoStream(object, lengthStream);
  return +lengthSum;
};

var coordinates = [null, null];
var object = {type: "LineString", coordinates: coordinates};

var distance = function(a, b) {
  coordinates[0] = a;
  coordinates[1] = b;
  return length(object);
};

var containsObjectType = {
  Feature: function(object, point) {
    return containsGeometry(object.geometry, point);
  },
  FeatureCollection: function(object, point) {
    var features = object.features, i = -1, n = features.length;
    while (++i < n) if (containsGeometry(features[i].geometry, point)) return true;
    return false;
  }
};

var containsGeometryType = {
  Sphere: function() {
    return true;
  },
  Point: function(object, point) {
    return containsPoint(object.coordinates, point);
  },
  MultiPoint: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPoint(coordinates[i], point)) return true;
    return false;
  },
  LineString: function(object, point) {
    return containsLine(object.coordinates, point);
  },
  MultiLineString: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsLine(coordinates[i], point)) return true;
    return false;
  },
  Polygon: function(object, point) {
    return containsPolygon(object.coordinates, point);
  },
  MultiPolygon: function(object, point) {
    var coordinates = object.coordinates, i = -1, n = coordinates.length;
    while (++i < n) if (containsPolygon(coordinates[i], point)) return true;
    return false;
  },
  GeometryCollection: function(object, point) {
    var geometries = object.geometries, i = -1, n = geometries.length;
    while (++i < n) if (containsGeometry(geometries[i], point)) return true;
    return false;
  }
};

function containsGeometry(geometry, point) {
  return geometry && containsGeometryType.hasOwnProperty(geometry.type)
      ? containsGeometryType[geometry.type](geometry, point)
      : false;
}

function containsPoint(coordinates, point) {
  return distance(coordinates, point) === 0;
}

function containsLine(coordinates, point) {
  var ab = distance(coordinates[0], coordinates[1]),
      ao = distance(coordinates[0], point),
      ob = distance(point, coordinates[1]);
  return ao + ob <= ab + epsilon;
}

function containsPolygon(coordinates, point) {
  return !!polygonContains(coordinates.map(ringRadians), pointRadians(point));
}

function ringRadians(ring) {
  return ring = ring.map(pointRadians), ring.pop(), ring;
}

function pointRadians(point) {
  return [point[0] * radians, point[1] * radians];
}

var contains = function(object, point) {
  return (object && containsObjectType.hasOwnProperty(object.type)
      ? containsObjectType[object.type]
      : containsGeometry)(object, point);
};

function graticuleX(y0, y1, dy) {
  var y = d3Array.range(y0, y1 - epsilon, dy).concat(y1);
  return function(x) { return y.map(function(y) { return [x, y]; }); };
}

function graticuleY(x0, x1, dx) {
  var x = d3Array.range(x0, x1 - epsilon, dx).concat(x1);
  return function(y) { return x.map(function(x) { return [x, y]; }); };
}

function graticule() {
  var x1, x0, X1, X0,
      y1, y0, Y1, Y0,
      dx = 10, dy = dx, DX = 90, DY = 360,
      x, y, X, Y,
      precision = 2.5;

  function graticule() {
    return {type: "MultiLineString", coordinates: lines()};
  }

  function lines() {
    return d3Array.range(ceil(X0 / DX) * DX, X1, DX).map(X)
        .concat(d3Array.range(ceil(Y0 / DY) * DY, Y1, DY).map(Y))
        .concat(d3Array.range(ceil(x0 / dx) * dx, x1, dx).filter(function(x) { return abs(x % DX) > epsilon; }).map(x))
        .concat(d3Array.range(ceil(y0 / dy) * dy, y1, dy).filter(function(y) { return abs(y % DY) > epsilon; }).map(y));
  }

  graticule.lines = function() {
    return lines().map(function(coordinates) { return {type: "LineString", coordinates: coordinates}; });
  };

  graticule.outline = function() {
    return {
      type: "Polygon",
      coordinates: [
        X(X0).concat(
        Y(Y1).slice(1),
        X(X1).reverse().slice(1),
        Y(Y0).reverse().slice(1))
      ]
    };
  };

  graticule.extent = function(_) {
    if (!arguments.length) return graticule.extentMinor();
    return graticule.extentMajor(_).extentMinor(_);
  };

  graticule.extentMajor = function(_) {
    if (!arguments.length) return [[X0, Y0], [X1, Y1]];
    X0 = +_[0][0], X1 = +_[1][0];
    Y0 = +_[0][1], Y1 = +_[1][1];
    if (X0 > X1) _ = X0, X0 = X1, X1 = _;
    if (Y0 > Y1) _ = Y0, Y0 = Y1, Y1 = _;
    return graticule.precision(precision);
  };

  graticule.extentMinor = function(_) {
    if (!arguments.length) return [[x0, y0], [x1, y1]];
    x0 = +_[0][0], x1 = +_[1][0];
    y0 = +_[0][1], y1 = +_[1][1];
    if (x0 > x1) _ = x0, x0 = x1, x1 = _;
    if (y0 > y1) _ = y0, y0 = y1, y1 = _;
    return graticule.precision(precision);
  };

  graticule.step = function(_) {
    if (!arguments.length) return graticule.stepMinor();
    return graticule.stepMajor(_).stepMinor(_);
  };

  graticule.stepMajor = function(_) {
    if (!arguments.length) return [DX, DY];
    DX = +_[0], DY = +_[1];
    return graticule;
  };

  graticule.stepMinor = function(_) {
    if (!arguments.length) return [dx, dy];
    dx = +_[0], dy = +_[1];
    return graticule;
  };

  graticule.precision = function(_) {
    if (!arguments.length) return precision;
    precision = +_;
    x = graticuleX(y0, y1, 90);
    y = graticuleY(x0, x1, precision);
    X = graticuleX(Y0, Y1, 90);
    Y = graticuleY(X0, X1, precision);
    return graticule;
  };

  return graticule
      .extentMajor([[-180, -90 + epsilon], [180, 90 - epsilon]])
      .extentMinor([[-180, -80 - epsilon], [180, 80 + epsilon]]);
}

function graticule10() {
  return graticule()();
}

var interpolate = function(a, b) {
  var x0 = a[0] * radians,
      y0 = a[1] * radians,
      x1 = b[0] * radians,
      y1 = b[1] * radians,
      cy0 = cos(y0),
      sy0 = sin(y0),
      cy1 = cos(y1),
      sy1 = sin(y1),
      kx0 = cy0 * cos(x0),
      ky0 = cy0 * sin(x0),
      kx1 = cy1 * cos(x1),
      ky1 = cy1 * sin(x1),
      d = 2 * asin(sqrt(haversin(y1 - y0) + cy0 * cy1 * haversin(x1 - x0))),
      k = sin(d);

  var interpolate = d ? function(t) {
    var B = sin(t *= d) / k,
        A = sin(d - t) / k,
        x = A * kx0 + B * kx1,
        y = A * ky0 + B * ky1,
        z = A * sy0 + B * sy1;
    return [
      atan2(y, x) * degrees,
      atan2(z, sqrt(x * x + y * y)) * degrees
    ];
  } : function() {
    return [x0 * degrees, y0 * degrees];
  };

  interpolate.distance = d;

  return interpolate;
};

var identity = function(x) {
  return x;
};

var areaSum$1 = adder();
var areaRingSum$1 = adder();
var x00;
var y00;
var x0$1;
var y0$1;

var areaStream$1 = {
  point: noop,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: function() {
    areaStream$1.lineStart = areaRingStart$1;
    areaStream$1.lineEnd = areaRingEnd$1;
  },
  polygonEnd: function() {
    areaStream$1.lineStart = areaStream$1.lineEnd = areaStream$1.point = noop;
    areaSum$1.add(abs(areaRingSum$1));
    areaRingSum$1.reset();
  },
  result: function() {
    var area = areaSum$1 / 2;
    areaSum$1.reset();
    return area;
  }
};

function areaRingStart$1() {
  areaStream$1.point = areaPointFirst$1;
}

function areaPointFirst$1(x, y) {
  areaStream$1.point = areaPoint$1;
  x00 = x0$1 = x, y00 = y0$1 = y;
}

function areaPoint$1(x, y) {
  areaRingSum$1.add(y0$1 * x - x0$1 * y);
  x0$1 = x, y0$1 = y;
}

function areaRingEnd$1() {
  areaPoint$1(x00, y00);
}

var x0$2 = Infinity;
var y0$2 = x0$2;
var x1 = -x0$2;
var y1 = x1;

var boundsStream$1 = {
  point: boundsPoint$1,
  lineStart: noop,
  lineEnd: noop,
  polygonStart: noop,
  polygonEnd: noop,
  result: function() {
    var bounds = [[x0$2, y0$2], [x1, y1]];
    x1 = y1 = -(y0$2 = x0$2 = Infinity);
    return bounds;
  }
};

function boundsPoint$1(x, y) {
  if (x < x0$2) x0$2 = x;
  if (x > x1) x1 = x;
  if (y < y0$2) y0$2 = y;
  if (y > y1) y1 = y;
}

// TODO Enforce positive area for exterior, negative area for interior?

var X0$1 = 0;
var Y0$1 = 0;
var Z0$1 = 0;
var X1$1 = 0;
var Y1$1 = 0;
var Z1$1 = 0;
var X2$1 = 0;
var Y2$1 = 0;
var Z2$1 = 0;
var x00$1;
var y00$1;
var x0$3;
var y0$3;

var centroidStream$1 = {
  point: centroidPoint$1,
  lineStart: centroidLineStart$1,
  lineEnd: centroidLineEnd$1,
  polygonStart: function() {
    centroidStream$1.lineStart = centroidRingStart$1;
    centroidStream$1.lineEnd = centroidRingEnd$1;
  },
  polygonEnd: function() {
    centroidStream$1.point = centroidPoint$1;
    centroidStream$1.lineStart = centroidLineStart$1;
    centroidStream$1.lineEnd = centroidLineEnd$1;
  },
  result: function() {
    var centroid = Z2$1 ? [X2$1 / Z2$1, Y2$1 / Z2$1]
        : Z1$1 ? [X1$1 / Z1$1, Y1$1 / Z1$1]
        : Z0$1 ? [X0$1 / Z0$1, Y0$1 / Z0$1]
        : [NaN, NaN];
    X0$1 = Y0$1 = Z0$1 =
    X1$1 = Y1$1 = Z1$1 =
    X2$1 = Y2$1 = Z2$1 = 0;
    return centroid;
  }
};

function centroidPoint$1(x, y) {
  X0$1 += x;
  Y0$1 += y;
  ++Z0$1;
}

function centroidLineStart$1() {
  centroidStream$1.point = centroidPointFirstLine;
}

function centroidPointFirstLine(x, y) {
  centroidStream$1.point = centroidPointLine;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function centroidPointLine(x, y) {
  var dx = x - x0$3, dy = y - y0$3, z = sqrt(dx * dx + dy * dy);
  X1$1 += z * (x0$3 + x) / 2;
  Y1$1 += z * (y0$3 + y) / 2;
  Z1$1 += z;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function centroidLineEnd$1() {
  centroidStream$1.point = centroidPoint$1;
}

function centroidRingStart$1() {
  centroidStream$1.point = centroidPointFirstRing;
}

function centroidRingEnd$1() {
  centroidPointRing(x00$1, y00$1);
}

function centroidPointFirstRing(x, y) {
  centroidStream$1.point = centroidPointRing;
  centroidPoint$1(x00$1 = x0$3 = x, y00$1 = y0$3 = y);
}

function centroidPointRing(x, y) {
  var dx = x - x0$3,
      dy = y - y0$3,
      z = sqrt(dx * dx + dy * dy);

  X1$1 += z * (x0$3 + x) / 2;
  Y1$1 += z * (y0$3 + y) / 2;
  Z1$1 += z;

  z = y0$3 * x - x0$3 * y;
  X2$1 += z * (x0$3 + x);
  Y2$1 += z * (y0$3 + y);
  Z2$1 += z * 3;
  centroidPoint$1(x0$3 = x, y0$3 = y);
}

function PathContext(context) {
  this._context = context;
}

PathContext.prototype = {
  _radius: 4.5,
  pointRadius: function(_) {
    return this._radius = _, this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._context.closePath();
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._context.moveTo(x, y);
        this._point = 1;
        break;
      }
      case 1: {
        this._context.lineTo(x, y);
        break;
      }
      default: {
        this._context.moveTo(x + this._radius, y);
        this._context.arc(x, y, this._radius, 0, tau);
        break;
      }
    }
  },
  result: noop
};

var lengthSum$1 = adder();
var lengthRing;
var x00$2;
var y00$2;
var x0$4;
var y0$4;

var lengthStream$1 = {
  point: noop,
  lineStart: function() {
    lengthStream$1.point = lengthPointFirst$1;
  },
  lineEnd: function() {
    if (lengthRing) lengthPoint$1(x00$2, y00$2);
    lengthStream$1.point = noop;
  },
  polygonStart: function() {
    lengthRing = true;
  },
  polygonEnd: function() {
    lengthRing = null;
  },
  result: function() {
    var length = +lengthSum$1;
    lengthSum$1.reset();
    return length;
  }
};

function lengthPointFirst$1(x, y) {
  lengthStream$1.point = lengthPoint$1;
  x00$2 = x0$4 = x, y00$2 = y0$4 = y;
}

function lengthPoint$1(x, y) {
  x0$4 -= x, y0$4 -= y;
  lengthSum$1.add(sqrt(x0$4 * x0$4 + y0$4 * y0$4));
  x0$4 = x, y0$4 = y;
}

function PathString() {
  this._string = [];
}

PathString.prototype = {
  _radius: 4.5,
  _circle: circle$1(4.5),
  pointRadius: function(_) {
    if ((_ = +_) !== this._radius) this._radius = _, this._circle = null;
    return this;
  },
  polygonStart: function() {
    this._line = 0;
  },
  polygonEnd: function() {
    this._line = NaN;
  },
  lineStart: function() {
    this._point = 0;
  },
  lineEnd: function() {
    if (this._line === 0) this._string.push("Z");
    this._point = NaN;
  },
  point: function(x, y) {
    switch (this._point) {
      case 0: {
        this._string.push("M", x, ",", y);
        this._point = 1;
        break;
      }
      case 1: {
        this._string.push("L", x, ",", y);
        break;
      }
      default: {
        if (this._circle == null) this._circle = circle$1(this._radius);
        this._string.push("M", x, ",", y, this._circle);
        break;
      }
    }
  },
  result: function() {
    if (this._string.length) {
      var result = this._string.join("");
      this._string = [];
      return result;
    } else {
      return null;
    }
  }
};

function circle$1(radius) {
  return "m0," + radius
      + "a" + radius + "," + radius + " 0 1,1 0," + -2 * radius
      + "a" + radius + "," + radius + " 0 1,1 0," + 2 * radius
      + "z";
}

var index = function(projection, context) {
  var pointRadius = 4.5,
      projectionStream,
      contextStream;

  function path(object) {
    if (object) {
      if (typeof pointRadius === "function") contextStream.pointRadius(+pointRadius.apply(this, arguments));
      geoStream(object, projectionStream(contextStream));
    }
    return contextStream.result();
  }

  path.area = function(object) {
    geoStream(object, projectionStream(areaStream$1));
    return areaStream$1.result();
  };

  path.measure = function(object) {
    geoStream(object, projectionStream(lengthStream$1));
    return lengthStream$1.result();
  };

  path.bounds = function(object) {
    geoStream(object, projectionStream(boundsStream$1));
    return boundsStream$1.result();
  };

  path.centroid = function(object) {
    geoStream(object, projectionStream(centroidStream$1));
    return centroidStream$1.result();
  };

  path.projection = function(_) {
    return arguments.length ? (projectionStream = _ == null ? (projection = null, identity) : (projection = _).stream, path) : projection;
  };

  path.context = function(_) {
    if (!arguments.length) return context;
    contextStream = _ == null ? (context = null, new PathString) : new PathContext(context = _);
    if (typeof pointRadius !== "function") contextStream.pointRadius(pointRadius);
    return path;
  };

  path.pointRadius = function(_) {
    if (!arguments.length) return pointRadius;
    pointRadius = typeof _ === "function" ? _ : (contextStream.pointRadius(+_), +_);
    return path;
  };

  return path.projection(projection).context(context);
};

var clip = function(pointVisible, clipLine, interpolate, start) {
  return function(rotate, sink) {
    var line = clipLine(sink),
        rotatedStart = rotate.invert(start[0], start[1]),
        ringBuffer = clipBuffer(),
        ringSink = clipLine(ringBuffer),
        polygonStarted = false,
        polygon,
        segments,
        ring;

    var clip = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() {
        clip.point = pointRing;
        clip.lineStart = ringStart;
        clip.lineEnd = ringEnd;
        segments = [];
        polygon = [];
      },
      polygonEnd: function() {
        clip.point = point;
        clip.lineStart = lineStart;
        clip.lineEnd = lineEnd;
        segments = d3Array.merge(segments);
        var startInside = polygonContains(polygon, rotatedStart);
        if (segments.length) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          clipPolygon(segments, compareIntersection, startInside, interpolate, sink);
        } else if (startInside) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          interpolate(null, null, 1, sink);
          sink.lineEnd();
        }
        if (polygonStarted) sink.polygonEnd(), polygonStarted = false;
        segments = polygon = null;
      },
      sphere: function() {
        sink.polygonStart();
        sink.lineStart();
        interpolate(null, null, 1, sink);
        sink.lineEnd();
        sink.polygonEnd();
      }
    };

    function point(lambda, phi) {
      var point = rotate(lambda, phi);
      if (pointVisible(lambda = point[0], phi = point[1])) sink.point(lambda, phi);
    }

    function pointLine(lambda, phi) {
      var point = rotate(lambda, phi);
      line.point(point[0], point[1]);
    }

    function lineStart() {
      clip.point = pointLine;
      line.lineStart();
    }

    function lineEnd() {
      clip.point = point;
      line.lineEnd();
    }

    function pointRing(lambda, phi) {
      ring.push([lambda, phi]);
      var point = rotate(lambda, phi);
      ringSink.point(point[0], point[1]);
    }

    function ringStart() {
      ringSink.lineStart();
      ring = [];
    }

    function ringEnd() {
      pointRing(ring[0][0], ring[0][1]);
      ringSink.lineEnd();

      var clean = ringSink.clean(),
          ringSegments = ringBuffer.result(),
          i, n = ringSegments.length, m,
          segment,
          point;

      ring.pop();
      polygon.push(ring);
      ring = null;

      if (!n) return;

      // No intersections.
      if (clean & 1) {
        segment = ringSegments[0];
        if ((m = segment.length - 1) > 0) {
          if (!polygonStarted) sink.polygonStart(), polygonStarted = true;
          sink.lineStart();
          for (i = 0; i < m; ++i) sink.point((point = segment[i])[0], point[1]);
          sink.lineEnd();
        }
        return;
      }

      // Rejoin connected segments.
      // TODO reuse ringBuffer.rejoin()?
      if (n > 1 && clean & 2) ringSegments.push(ringSegments.pop().concat(ringSegments.shift()));

      segments.push(ringSegments.filter(validSegment));
    }

    return clip;
  };
};

function validSegment(segment) {
  return segment.length > 1;
}

// Intersections are sorted along the clip edge. For both antimeridian cutting
// and circle clipping, the same comparison is used.
function compareIntersection(a, b) {
  return ((a = a.x)[0] < 0 ? a[1] - halfPi - epsilon : halfPi - a[1])
       - ((b = b.x)[0] < 0 ? b[1] - halfPi - epsilon : halfPi - b[1]);
}

var clipAntimeridian = clip(
  function() { return true; },
  clipAntimeridianLine,
  clipAntimeridianInterpolate,
  [-pi, -halfPi]
);

// Takes a line and cuts into visible segments. Return values: 0 - there were
// intersections or the line was empty; 1 - no intersections; 2 - there were
// intersections, and the first and last segments should be rejoined.
function clipAntimeridianLine(stream) {
  var lambda0 = NaN,
      phi0 = NaN,
      sign0 = NaN,
      clean; // no intersections

  return {
    lineStart: function() {
      stream.lineStart();
      clean = 1;
    },
    point: function(lambda1, phi1) {
      var sign1 = lambda1 > 0 ? pi : -pi,
          delta = abs(lambda1 - lambda0);
      if (abs(delta - pi) < epsilon) { // line crosses a pole
        stream.point(lambda0, phi0 = (phi0 + phi1) / 2 > 0 ? halfPi : -halfPi);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        stream.point(lambda1, phi0);
        clean = 0;
      } else if (sign0 !== sign1 && delta >= pi) { // line crosses antimeridian
        if (abs(lambda0 - sign0) < epsilon) lambda0 -= sign0 * epsilon; // handle degeneracies
        if (abs(lambda1 - sign1) < epsilon) lambda1 -= sign1 * epsilon;
        phi0 = clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1);
        stream.point(sign0, phi0);
        stream.lineEnd();
        stream.lineStart();
        stream.point(sign1, phi0);
        clean = 0;
      }
      stream.point(lambda0 = lambda1, phi0 = phi1);
      sign0 = sign1;
    },
    lineEnd: function() {
      stream.lineEnd();
      lambda0 = phi0 = NaN;
    },
    clean: function() {
      return 2 - clean; // if intersections, rejoin first and last segments
    }
  };
}

function clipAntimeridianIntersect(lambda0, phi0, lambda1, phi1) {
  var cosPhi0,
      cosPhi1,
      sinLambda0Lambda1 = sin(lambda0 - lambda1);
  return abs(sinLambda0Lambda1) > epsilon
      ? atan((sin(phi0) * (cosPhi1 = cos(phi1)) * sin(lambda1)
          - sin(phi1) * (cosPhi0 = cos(phi0)) * sin(lambda0))
          / (cosPhi0 * cosPhi1 * sinLambda0Lambda1))
      : (phi0 + phi1) / 2;
}

function clipAntimeridianInterpolate(from, to, direction, stream) {
  var phi;
  if (from == null) {
    phi = direction * halfPi;
    stream.point(-pi, phi);
    stream.point(0, phi);
    stream.point(pi, phi);
    stream.point(pi, 0);
    stream.point(pi, -phi);
    stream.point(0, -phi);
    stream.point(-pi, -phi);
    stream.point(-pi, 0);
    stream.point(-pi, phi);
  } else if (abs(from[0] - to[0]) > epsilon) {
    var lambda = from[0] < to[0] ? pi : -pi;
    phi = direction * lambda / 2;
    stream.point(-lambda, phi);
    stream.point(0, phi);
    stream.point(lambda, phi);
  } else {
    stream.point(to[0], to[1]);
  }
}

var clipCircle = function(radius, delta) {
  var cr = cos(radius),
      smallRadius = cr > 0,
      notHemisphere = abs(cr) > epsilon; // TODO optimise for this common case

  function interpolate(from, to, direction, stream) {
    circleStream(stream, radius, delta, direction, from, to);
  }

  function visible(lambda, phi) {
    return cos(lambda) * cos(phi) > cr;
  }

  // Takes a line and cuts into visible segments. Return values used for polygon
  // clipping: 0 - there were intersections or the line was empty; 1 - no
  // intersections 2 - there were intersections, and the first and last segments
  // should be rejoined.
  function clipLine(stream) {
    var point0, // previous point
        c0, // code for previous point
        v0, // visibility of previous point
        v00, // visibility of first point
        clean; // no intersections
    return {
      lineStart: function() {
        v00 = v0 = false;
        clean = 1;
      },
      point: function(lambda, phi) {
        var point1 = [lambda, phi],
            point2,
            v = visible(lambda, phi),
            c = smallRadius
              ? v ? 0 : code(lambda, phi)
              : v ? code(lambda + (lambda < 0 ? pi : -pi), phi) : 0;
        if (!point0 && (v00 = v0 = v)) stream.lineStart();
        // Handle degeneracies.
        // TODO ignore if not clipping polygons.
        if (v !== v0) {
          point2 = intersect(point0, point1);
          if (!point2 || pointEqual(point0, point2) || pointEqual(point1, point2)) {
            point1[0] += epsilon;
            point1[1] += epsilon;
            v = visible(point1[0], point1[1]);
          }
        }
        if (v !== v0) {
          clean = 0;
          if (v) {
            // outside going in
            stream.lineStart();
            point2 = intersect(point1, point0);
            stream.point(point2[0], point2[1]);
          } else {
            // inside going out
            point2 = intersect(point0, point1);
            stream.point(point2[0], point2[1]);
            stream.lineEnd();
          }
          point0 = point2;
        } else if (notHemisphere && point0 && smallRadius ^ v) {
          var t;
          // If the codes for two points are different, or are both zero,
          // and there this segment intersects with the small circle.
          if (!(c & c0) && (t = intersect(point1, point0, true))) {
            clean = 0;
            if (smallRadius) {
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
            } else {
              stream.point(t[1][0], t[1][1]);
              stream.lineEnd();
              stream.lineStart();
              stream.point(t[0][0], t[0][1]);
            }
          }
        }
        if (v && (!point0 || !pointEqual(point0, point1))) {
          stream.point(point1[0], point1[1]);
        }
        point0 = point1, v0 = v, c0 = c;
      },
      lineEnd: function() {
        if (v0) stream.lineEnd();
        point0 = null;
      },
      // Rejoin first and last segments if there were intersections and the first
      // and last points were visible.
      clean: function() {
        return clean | ((v00 && v0) << 1);
      }
    };
  }

  // Intersects the great circle between a and b with the clip circle.
  function intersect(a, b, two) {
    var pa = cartesian(a),
        pb = cartesian(b);

    // We have two planes, n1.p = d1 and n2.p = d2.
    // Find intersection line p(t) = c1 n1 + c2 n2 + t (n1 ⨯ n2).
    var n1 = [1, 0, 0], // normal
        n2 = cartesianCross(pa, pb),
        n2n2 = cartesianDot(n2, n2),
        n1n2 = n2[0], // cartesianDot(n1, n2),
        determinant = n2n2 - n1n2 * n1n2;

    // Two polar points.
    if (!determinant) return !two && a;

    var c1 =  cr * n2n2 / determinant,
        c2 = -cr * n1n2 / determinant,
        n1xn2 = cartesianCross(n1, n2),
        A = cartesianScale(n1, c1),
        B = cartesianScale(n2, c2);
    cartesianAddInPlace(A, B);

    // Solve |p(t)|^2 = 1.
    var u = n1xn2,
        w = cartesianDot(A, u),
        uu = cartesianDot(u, u),
        t2 = w * w - uu * (cartesianDot(A, A) - 1);

    if (t2 < 0) return;

    var t = sqrt(t2),
        q = cartesianScale(u, (-w - t) / uu);
    cartesianAddInPlace(q, A);
    q = spherical(q);

    if (!two) return q;

    // Two intersection points.
    var lambda0 = a[0],
        lambda1 = b[0],
        phi0 = a[1],
        phi1 = b[1],
        z;

    if (lambda1 < lambda0) z = lambda0, lambda0 = lambda1, lambda1 = z;

    var delta = lambda1 - lambda0,
        polar = abs(delta - pi) < epsilon,
        meridian = polar || delta < epsilon;

    if (!polar && phi1 < phi0) z = phi0, phi0 = phi1, phi1 = z;

    // Check that the first point is between a and b.
    if (meridian
        ? polar
          ? phi0 + phi1 > 0 ^ q[1] < (abs(q[0] - lambda0) < epsilon ? phi0 : phi1)
          : phi0 <= q[1] && q[1] <= phi1
        : delta > pi ^ (lambda0 <= q[0] && q[0] <= lambda1)) {
      var q1 = cartesianScale(u, (-w + t) / uu);
      cartesianAddInPlace(q1, A);
      return [q, spherical(q1)];
    }
  }

  // Generates a 4-bit vector representing the location of a point relative to
  // the small circle's bounding box.
  function code(lambda, phi) {
    var r = smallRadius ? radius : pi - radius,
        code = 0;
    if (lambda < -r) code |= 1; // left
    else if (lambda > r) code |= 2; // right
    if (phi < -r) code |= 4; // below
    else if (phi > r) code |= 8; // above
    return code;
  }

  return clip(visible, clipLine, interpolate, smallRadius ? [0, -radius] : [-pi, radius - pi]);
};

var transform = function(methods) {
  return {
    stream: transformer(methods)
  };
};

function transformer(methods) {
  return function(stream) {
    var s = new TransformStream;
    for (var key in methods) s[key] = methods[key];
    s.stream = stream;
    return s;
  };
}

function TransformStream() {}

TransformStream.prototype = {
  constructor: TransformStream,
  point: function(x, y) { this.stream.point(x, y); },
  sphere: function() { this.stream.sphere(); },
  lineStart: function() { this.stream.lineStart(); },
  lineEnd: function() { this.stream.lineEnd(); },
  polygonStart: function() { this.stream.polygonStart(); },
  polygonEnd: function() { this.stream.polygonEnd(); }
};

function fitExtent(projection, extent, object) {
  var w = extent[1][0] - extent[0][0],
      h = extent[1][1] - extent[0][1],
      clip = projection.clipExtent && projection.clipExtent();

  projection
      .scale(150)
      .translate([0, 0]);

  if (clip != null) projection.clipExtent(null);

  geoStream(object, projection.stream(boundsStream$1));

  var b = boundsStream$1.result(),
      k = Math.min(w / (b[1][0] - b[0][0]), h / (b[1][1] - b[0][1])),
      x = +extent[0][0] + (w - k * (b[1][0] + b[0][0])) / 2,
      y = +extent[0][1] + (h - k * (b[1][1] + b[0][1])) / 2;

  if (clip != null) projection.clipExtent(clip);

  return projection
      .scale(k * 150)
      .translate([x, y]);
}

function fitSize(projection, size, object) {
  return fitExtent(projection, [[0, 0], size], object);
}

var maxDepth = 16;
var cosMinDistance = cos(30 * radians); // cos(minimum angular distance)

var resample = function(project, delta2) {
  return +delta2 ? resample$1(project, delta2) : resampleNone(project);
};

function resampleNone(project) {
  return transformer({
    point: function(x, y) {
      x = project(x, y);
      this.stream.point(x[0], x[1]);
    }
  });
}

function resample$1(project, delta2) {

  function resampleLineTo(x0, y0, lambda0, a0, b0, c0, x1, y1, lambda1, a1, b1, c1, depth, stream) {
    var dx = x1 - x0,
        dy = y1 - y0,
        d2 = dx * dx + dy * dy;
    if (d2 > 4 * delta2 && depth--) {
      var a = a0 + a1,
          b = b0 + b1,
          c = c0 + c1,
          m = sqrt(a * a + b * b + c * c),
          phi2 = asin(c /= m),
          lambda2 = abs(abs(c) - 1) < epsilon || abs(lambda0 - lambda1) < epsilon ? (lambda0 + lambda1) / 2 : atan2(b, a),
          p = project(lambda2, phi2),
          x2 = p[0],
          y2 = p[1],
          dx2 = x2 - x0,
          dy2 = y2 - y0,
          dz = dy * dx2 - dx * dy2;
      if (dz * dz / d2 > delta2 // perpendicular projected distance
          || abs((dx * dx2 + dy * dy2) / d2 - 0.5) > 0.3 // midpoint close to an end
          || a0 * a1 + b0 * b1 + c0 * c1 < cosMinDistance) { // angular distance
        resampleLineTo(x0, y0, lambda0, a0, b0, c0, x2, y2, lambda2, a /= m, b /= m, c, depth, stream);
        stream.point(x2, y2);
        resampleLineTo(x2, y2, lambda2, a, b, c, x1, y1, lambda1, a1, b1, c1, depth, stream);
      }
    }
  }
  return function(stream) {
    var lambda00, x00, y00, a00, b00, c00, // first point
        lambda0, x0, y0, a0, b0, c0; // previous point

    var resampleStream = {
      point: point,
      lineStart: lineStart,
      lineEnd: lineEnd,
      polygonStart: function() { stream.polygonStart(); resampleStream.lineStart = ringStart; },
      polygonEnd: function() { stream.polygonEnd(); resampleStream.lineStart = lineStart; }
    };

    function point(x, y) {
      x = project(x, y);
      stream.point(x[0], x[1]);
    }

    function lineStart() {
      x0 = NaN;
      resampleStream.point = linePoint;
      stream.lineStart();
    }

    function linePoint(lambda, phi) {
      var c = cartesian([lambda, phi]), p = project(lambda, phi);
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x0 = p[0], y0 = p[1], lambda0 = lambda, a0 = c[0], b0 = c[1], c0 = c[2], maxDepth, stream);
      stream.point(x0, y0);
    }

    function lineEnd() {
      resampleStream.point = point;
      stream.lineEnd();
    }

    function ringStart() {
      lineStart();
      resampleStream.point = ringPoint;
      resampleStream.lineEnd = ringEnd;
    }

    function ringPoint(lambda, phi) {
      linePoint(lambda00 = lambda, phi), x00 = x0, y00 = y0, a00 = a0, b00 = b0, c00 = c0;
      resampleStream.point = linePoint;
    }

    function ringEnd() {
      resampleLineTo(x0, y0, lambda0, a0, b0, c0, x00, y00, lambda00, a00, b00, c00, maxDepth, stream);
      resampleStream.lineEnd = lineEnd;
      lineEnd();
    }

    return resampleStream;
  };
}

var transformRadians = transformer({
  point: function(x, y) {
    this.stream.point(x * radians, y * radians);
  }
});

function projection(project) {
  return projectionMutator(function() { return project; })();
}

function projectionMutator(projectAt) {
  var project,
      k = 150, // scale
      x = 480, y = 250, // translate
      dx, dy, lambda = 0, phi = 0, // center
      deltaLambda = 0, deltaPhi = 0, deltaGamma = 0, rotate, projectRotate, // rotate
      theta = null, preclip = clipAntimeridian, // clip angle
      x0 = null, y0, x1, y1, postclip = identity, // clip extent
      delta2 = 0.5, projectResample = resample(projectTransform, delta2), // precision
      cache,
      cacheStream;

  function projection(point) {
    point = projectRotate(point[0] * radians, point[1] * radians);
    return [point[0] * k + dx, dy - point[1] * k];
  }

  function invert(point) {
    point = projectRotate.invert((point[0] - dx) / k, (dy - point[1]) / k);
    return point && [point[0] * degrees, point[1] * degrees];
  }

  function projectTransform(x, y) {
    return x = project(x, y), [x[0] * k + dx, dy - x[1] * k];
  }

  projection.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = transformRadians(preclip(rotate, projectResample(postclip(cacheStream = stream))));
  };

  projection.clipAngle = function(_) {
    return arguments.length ? (preclip = +_ ? clipCircle(theta = _ * radians, 6 * radians) : (theta = null, clipAntimeridian), reset()) : theta * degrees;
  };

  projection.clipExtent = function(_) {
    return arguments.length ? (postclip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  projection.scale = function(_) {
    return arguments.length ? (k = +_, recenter()) : k;
  };

  projection.translate = function(_) {
    return arguments.length ? (x = +_[0], y = +_[1], recenter()) : [x, y];
  };

  projection.center = function(_) {
    return arguments.length ? (lambda = _[0] % 360 * radians, phi = _[1] % 360 * radians, recenter()) : [lambda * degrees, phi * degrees];
  };

  projection.rotate = function(_) {
    return arguments.length ? (deltaLambda = _[0] % 360 * radians, deltaPhi = _[1] % 360 * radians, deltaGamma = _.length > 2 ? _[2] % 360 * radians : 0, recenter()) : [deltaLambda * degrees, deltaPhi * degrees, deltaGamma * degrees];
  };

  projection.precision = function(_) {
    return arguments.length ? (projectResample = resample(projectTransform, delta2 = _ * _), reset()) : sqrt(delta2);
  };

  projection.fitExtent = function(extent$$1, object) {
    return fitExtent(projection, extent$$1, object);
  };

  projection.fitSize = function(size, object) {
    return fitSize(projection, size, object);
  };

  function recenter() {
    projectRotate = compose(rotate = rotateRadians(deltaLambda, deltaPhi, deltaGamma), project);
    var center = project(lambda, phi);
    dx = x - center[0] * k;
    dy = y + center[1] * k;
    return reset();
  }

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return function() {
    project = projectAt.apply(this, arguments);
    projection.invert = project.invert && invert;
    return recenter();
  };
}

function conicProjection(projectAt) {
  var phi0 = 0,
      phi1 = pi / 3,
      m = projectionMutator(projectAt),
      p = m(phi0, phi1);

  p.parallels = function(_) {
    return arguments.length ? m(phi0 = _[0] * radians, phi1 = _[1] * radians) : [phi0 * degrees, phi1 * degrees];
  };

  return p;
}

function cylindricalEqualAreaRaw(phi0) {
  var cosPhi0 = cos(phi0);

  function forward(lambda, phi) {
    return [lambda * cosPhi0, sin(phi) / cosPhi0];
  }

  forward.invert = function(x, y) {
    return [x / cosPhi0, asin(y * cosPhi0)];
  };

  return forward;
}

function conicEqualAreaRaw(y0, y1) {
  var sy0 = sin(y0), n = (sy0 + sin(y1)) / 2;

  // Are the parallels symmetrical around the Equator?
  if (abs(n) < epsilon) return cylindricalEqualAreaRaw(y0);

  var c = 1 + sy0 * (2 * n - sy0), r0 = sqrt(c) / n;

  function project(x, y) {
    var r = sqrt(c - 2 * n * sin(y)) / n;
    return [r * sin(x *= n), r0 - r * cos(x)];
  }

  project.invert = function(x, y) {
    var r0y = r0 - y;
    return [atan2(x, abs(r0y)) / n * sign(r0y), asin((c - (x * x + r0y * r0y) * n * n) / (2 * n))];
  };

  return project;
}

var conicEqualArea = function() {
  return conicProjection(conicEqualAreaRaw)
      .scale(155.424)
      .center([0, 33.6442]);
};

var albers = function() {
  return conicEqualArea()
      .parallels([29.5, 45.5])
      .scale(1070)
      .translate([480, 250])
      .rotate([96, 0])
      .center([-0.6, 38.7]);
};

// The projections must have mutually exclusive clip regions on the sphere,
// as this will avoid emitting interleaving lines and polygons.
function multiplex(streams) {
  var n = streams.length;
  return {
    point: function(x, y) { var i = -1; while (++i < n) streams[i].point(x, y); },
    sphere: function() { var i = -1; while (++i < n) streams[i].sphere(); },
    lineStart: function() { var i = -1; while (++i < n) streams[i].lineStart(); },
    lineEnd: function() { var i = -1; while (++i < n) streams[i].lineEnd(); },
    polygonStart: function() { var i = -1; while (++i < n) streams[i].polygonStart(); },
    polygonEnd: function() { var i = -1; while (++i < n) streams[i].polygonEnd(); }
  };
}

// A composite projection for the United States, configured by default for
// 960×500. The projection also works quite well at 960×600 if you change the
// scale to 1285 and adjust the translate accordingly. The set of standard
// parallels for each region comes from USGS, which is published here:
// http://egsc.usgs.gov/isb/pubs/MapProjections/projections.html#albers
var albersUsa = function() {
  var cache,
      cacheStream,
      lower48 = albers(), lower48Point,
      alaska = conicEqualArea().rotate([154, 0]).center([-2, 58.5]).parallels([55, 65]), alaskaPoint, // EPSG:3338
      hawaii = conicEqualArea().rotate([157, 0]).center([-3, 19.9]).parallels([8, 18]), hawaiiPoint, // ESRI:102007
      point, pointStream = {point: function(x, y) { point = [x, y]; }};

  function albersUsa(coordinates) {
    var x = coordinates[0], y = coordinates[1];
    return point = null,
        (lower48Point.point(x, y), point)
        || (alaskaPoint.point(x, y), point)
        || (hawaiiPoint.point(x, y), point);
  }

  albersUsa.invert = function(coordinates) {
    var k = lower48.scale(),
        t = lower48.translate(),
        x = (coordinates[0] - t[0]) / k,
        y = (coordinates[1] - t[1]) / k;
    return (y >= 0.120 && y < 0.234 && x >= -0.425 && x < -0.214 ? alaska
        : y >= 0.166 && y < 0.234 && x >= -0.214 && x < -0.115 ? hawaii
        : lower48).invert(coordinates);
  };

  albersUsa.stream = function(stream) {
    return cache && cacheStream === stream ? cache : cache = multiplex([lower48.stream(cacheStream = stream), alaska.stream(stream), hawaii.stream(stream)]);
  };

  albersUsa.precision = function(_) {
    if (!arguments.length) return lower48.precision();
    lower48.precision(_), alaska.precision(_), hawaii.precision(_);
    return reset();
  };

  albersUsa.scale = function(_) {
    if (!arguments.length) return lower48.scale();
    lower48.scale(_), alaska.scale(_ * 0.35), hawaii.scale(_);
    return albersUsa.translate(lower48.translate());
  };

  albersUsa.translate = function(_) {
    if (!arguments.length) return lower48.translate();
    var k = lower48.scale(), x = +_[0], y = +_[1];

    lower48Point = lower48
        .translate(_)
        .clipExtent([[x - 0.455 * k, y - 0.238 * k], [x + 0.455 * k, y + 0.238 * k]])
        .stream(pointStream);

    alaskaPoint = alaska
        .translate([x - 0.307 * k, y + 0.201 * k])
        .clipExtent([[x - 0.425 * k + epsilon, y + 0.120 * k + epsilon], [x - 0.214 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    hawaiiPoint = hawaii
        .translate([x - 0.205 * k, y + 0.212 * k])
        .clipExtent([[x - 0.214 * k + epsilon, y + 0.166 * k + epsilon], [x - 0.115 * k - epsilon, y + 0.234 * k - epsilon]])
        .stream(pointStream);

    return reset();
  };

  albersUsa.fitExtent = function(extent, object) {
    return fitExtent(albersUsa, extent, object);
  };

  albersUsa.fitSize = function(size, object) {
    return fitSize(albersUsa, size, object);
  };

  function reset() {
    cache = cacheStream = null;
    return albersUsa;
  }

  return albersUsa.scale(1070);
};

function azimuthalRaw(scale) {
  return function(x, y) {
    var cx = cos(x),
        cy = cos(y),
        k = scale(cx * cy);
    return [
      k * cy * sin(x),
      k * sin(y)
    ];
  }
}

function azimuthalInvert(angle) {
  return function(x, y) {
    var z = sqrt(x * x + y * y),
        c = angle(z),
        sc = sin(c),
        cc = cos(c);
    return [
      atan2(x * sc, z * cc),
      asin(z && y * sc / z)
    ];
  }
}

var azimuthalEqualAreaRaw = azimuthalRaw(function(cxcy) {
  return sqrt(2 / (1 + cxcy));
});

azimuthalEqualAreaRaw.invert = azimuthalInvert(function(z) {
  return 2 * asin(z / 2);
});

var azimuthalEqualArea = function() {
  return projection(azimuthalEqualAreaRaw)
      .scale(124.75)
      .clipAngle(180 - 1e-3);
};

var azimuthalEquidistantRaw = azimuthalRaw(function(c) {
  return (c = acos(c)) && c / sin(c);
});

azimuthalEquidistantRaw.invert = azimuthalInvert(function(z) {
  return z;
});

var azimuthalEquidistant = function() {
  return projection(azimuthalEquidistantRaw)
      .scale(79.4188)
      .clipAngle(180 - 1e-3);
};

function mercatorRaw(lambda, phi) {
  return [lambda, log(tan((halfPi + phi) / 2))];
}

mercatorRaw.invert = function(x, y) {
  return [x, 2 * atan(exp(y)) - halfPi];
};

var mercator = function() {
  return mercatorProjection(mercatorRaw)
      .scale(961 / tau);
};

function mercatorProjection(project) {
  var m = projection(project),
      center = m.center,
      scale = m.scale,
      translate = m.translate,
      clipExtent = m.clipExtent,
      x0 = null, y0, x1, y1; // clip extent

  m.scale = function(_) {
    return arguments.length ? (scale(_), reclip()) : scale();
  };

  m.translate = function(_) {
    return arguments.length ? (translate(_), reclip()) : translate();
  };

  m.center = function(_) {
    return arguments.length ? (center(_), reclip()) : center();
  };

  m.clipExtent = function(_) {
    return arguments.length ? ((_ == null ? x0 = y0 = x1 = y1 = null : (x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1])), reclip()) : x0 == null ? null : [[x0, y0], [x1, y1]];
  };

  function reclip() {
    var k = pi * scale(),
        t = m(rotation(m.rotate()).invert([0, 0]));
    return clipExtent(x0 == null
        ? [[t[0] - k, t[1] - k], [t[0] + k, t[1] + k]] : project === mercatorRaw
        ? [[Math.max(t[0] - k, x0), y0], [Math.min(t[0] + k, x1), y1]]
        : [[x0, Math.max(t[1] - k, y0)], [x1, Math.min(t[1] + k, y1)]]);
  }

  return reclip();
}

function tany(y) {
  return tan((halfPi + y) / 2);
}

function conicConformalRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : log(cy0 / cos(y1)) / log(tany(y1) / tany(y0)),
      f = cy0 * pow(tany(y0), n) / n;

  if (!n) return mercatorRaw;

  function project(x, y) {
    if (f > 0) { if (y < -halfPi + epsilon) y = -halfPi + epsilon; }
    else { if (y > halfPi - epsilon) y = halfPi - epsilon; }
    var r = f / pow(tany(y), n);
    return [r * sin(n * x), f - r * cos(n * x)];
  }

  project.invert = function(x, y) {
    var fy = f - y, r = sign(n) * sqrt(x * x + fy * fy);
    return [atan2(x, abs(fy)) / n * sign(fy), 2 * atan(pow(f / r, 1 / n)) - halfPi];
  };

  return project;
}

var conicConformal = function() {
  return conicProjection(conicConformalRaw)
      .scale(109.5)
      .parallels([30, 30]);
};

function equirectangularRaw(lambda, phi) {
  return [lambda, phi];
}

equirectangularRaw.invert = equirectangularRaw;

var equirectangular = function() {
  return projection(equirectangularRaw)
      .scale(152.63);
};

function conicEquidistantRaw(y0, y1) {
  var cy0 = cos(y0),
      n = y0 === y1 ? sin(y0) : (cy0 - cos(y1)) / (y1 - y0),
      g = cy0 / n + y0;

  if (abs(n) < epsilon) return equirectangularRaw;

  function project(x, y) {
    var gy = g - y, nx = n * x;
    return [gy * sin(nx), g - gy * cos(nx)];
  }

  project.invert = function(x, y) {
    var gy = g - y;
    return [atan2(x, abs(gy)) / n * sign(gy), g - sign(n) * sqrt(x * x + gy * gy)];
  };

  return project;
}

var conicEquidistant = function() {
  return conicProjection(conicEquidistantRaw)
      .scale(131.154)
      .center([0, 13.9389]);
};

function gnomonicRaw(x, y) {
  var cy = cos(y), k = cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

gnomonicRaw.invert = azimuthalInvert(atan);

var gnomonic = function() {
  return projection(gnomonicRaw)
      .scale(144.049)
      .clipAngle(60);
};

function scaleTranslate(kx, ky, tx, ty) {
  return kx === 1 && ky === 1 && tx === 0 && ty === 0 ? identity : transformer({
    point: function(x, y) {
      this.stream.point(x * kx + tx, y * ky + ty);
    }
  });
}

var identity$1 = function() {
  var k = 1, tx = 0, ty = 0, sx = 1, sy = 1, transform$$1 = identity, // scale, translate and reflect
      x0 = null, y0, x1, y1, clip = identity, // clip extent
      cache,
      cacheStream,
      projection;

  function reset() {
    cache = cacheStream = null;
    return projection;
  }

  return projection = {
    stream: function(stream) {
      return cache && cacheStream === stream ? cache : cache = transform$$1(clip(cacheStream = stream));
    },
    clipExtent: function(_) {
      return arguments.length ? (clip = _ == null ? (x0 = y0 = x1 = y1 = null, identity) : clipExtent(x0 = +_[0][0], y0 = +_[0][1], x1 = +_[1][0], y1 = +_[1][1]), reset()) : x0 == null ? null : [[x0, y0], [x1, y1]];
    },
    scale: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate((k = +_) * sx, k * sy, tx, ty), reset()) : k;
    },
    translate: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate(k * sx, k * sy, tx = +_[0], ty = +_[1]), reset()) : [tx, ty];
    },
    reflectX: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate(k * (sx = _ ? -1 : 1), k * sy, tx, ty), reset()) : sx < 0;
    },
    reflectY: function(_) {
      return arguments.length ? (transform$$1 = scaleTranslate(k * sx, k * (sy = _ ? -1 : 1), tx, ty), reset()) : sy < 0;
    },
    fitExtent: function(extent$$1, object) {
      return fitExtent(projection, extent$$1, object);
    },
    fitSize: function(size, object) {
      return fitSize(projection, size, object);
    }
  };
};

function naturalEarth1Raw(lambda, phi) {
  var phi2 = phi * phi, phi4 = phi2 * phi2;
  return [
    lambda * (0.8707 - 0.131979 * phi2 + phi4 * (-0.013791 + phi4 * (0.003971 * phi2 - 0.001529 * phi4))),
    phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4)))
  ];
}

naturalEarth1Raw.invert = function(x, y) {
  var phi = y, i = 25, delta;
  do {
    var phi2 = phi * phi, phi4 = phi2 * phi2;
    phi -= delta = (phi * (1.007226 + phi2 * (0.015085 + phi4 * (-0.044475 + 0.028874 * phi2 - 0.005916 * phi4))) - y) /
        (1.007226 + phi2 * (0.015085 * 3 + phi4 * (-0.044475 * 7 + 0.028874 * 9 * phi2 - 0.005916 * 11 * phi4)));
  } while (abs(delta) > epsilon && --i > 0);
  return [
    x / (0.8707 + (phi2 = phi * phi) * (-0.131979 + phi2 * (-0.013791 + phi2 * phi2 * phi2 * (0.003971 - 0.001529 * phi2)))),
    phi
  ];
};

var naturalEarth1 = function() {
  return projection(naturalEarth1Raw)
      .scale(175.295);
};

function orthographicRaw(x, y) {
  return [cos(y) * sin(x), sin(y)];
}

orthographicRaw.invert = azimuthalInvert(asin);

var orthographic = function() {
  return projection(orthographicRaw)
      .scale(249.5)
      .clipAngle(90 + epsilon);
};

function stereographicRaw(x, y) {
  var cy = cos(y), k = 1 + cos(x) * cy;
  return [cy * sin(x) / k, sin(y) / k];
}

stereographicRaw.invert = azimuthalInvert(function(z) {
  return 2 * atan(z);
});

var stereographic = function() {
  return projection(stereographicRaw)
      .scale(250)
      .clipAngle(142);
};

function transverseMercatorRaw(lambda, phi) {
  return [log(tan((halfPi + phi) / 2)), -lambda];
}

transverseMercatorRaw.invert = function(x, y) {
  return [-y, 2 * atan(exp(x)) - halfPi];
};

var transverseMercator = function() {
  var m = mercatorProjection(transverseMercatorRaw),
      center = m.center,
      rotate = m.rotate;

  m.center = function(_) {
    return arguments.length ? center([-_[1], _[0]]) : (_ = center(), [_[1], -_[0]]);
  };

  m.rotate = function(_) {
    return arguments.length ? rotate([_[0], _[1], _.length > 2 ? _[2] + 90 : 90]) : (_ = rotate(), [_[0], _[1], _[2] - 90]);
  };

  return rotate([0, 0, 90])
      .scale(159.155);
};

exports.geoArea = area;
exports.geoBounds = bounds;
exports.geoCentroid = centroid;
exports.geoCircle = circle;
exports.geoClipExtent = extent;
exports.geoContains = contains;
exports.geoDistance = distance;
exports.geoGraticule = graticule;
exports.geoGraticule10 = graticule10;
exports.geoInterpolate = interpolate;
exports.geoLength = length;
exports.geoPath = index;
exports.geoAlbers = albers;
exports.geoAlbersUsa = albersUsa;
exports.geoAzimuthalEqualArea = azimuthalEqualArea;
exports.geoAzimuthalEqualAreaRaw = azimuthalEqualAreaRaw;
exports.geoAzimuthalEquidistant = azimuthalEquidistant;
exports.geoAzimuthalEquidistantRaw = azimuthalEquidistantRaw;
exports.geoConicConformal = conicConformal;
exports.geoConicConformalRaw = conicConformalRaw;
exports.geoConicEqualArea = conicEqualArea;
exports.geoConicEqualAreaRaw = conicEqualAreaRaw;
exports.geoConicEquidistant = conicEquidistant;
exports.geoConicEquidistantRaw = conicEquidistantRaw;
exports.geoEquirectangular = equirectangular;
exports.geoEquirectangularRaw = equirectangularRaw;
exports.geoGnomonic = gnomonic;
exports.geoGnomonicRaw = gnomonicRaw;
exports.geoIdentity = identity$1;
exports.geoProjection = projection;
exports.geoProjectionMutator = projectionMutator;
exports.geoMercator = mercator;
exports.geoMercatorRaw = mercatorRaw;
exports.geoNaturalEarth1 = naturalEarth1;
exports.geoNaturalEarth1Raw = naturalEarth1Raw;
exports.geoOrthographic = orthographic;
exports.geoOrthographicRaw = orthographicRaw;
exports.geoStereographic = stereographic;
exports.geoStereographicRaw = stereographicRaw;
exports.geoTransverseMercator = transverseMercator;
exports.geoTransverseMercatorRaw = transverseMercatorRaw;
exports.geoRotation = rotation;
exports.geoStream = geoStream;
exports.geoTransform = transform;

Object.defineProperty(exports, '__esModule', { value: true });

})));

},{"d3-array":49}],51:[function(require,module,exports){
var rbush = require('rbush');
var helpers = require('@turf/helpers');
var meta = require('@turf/meta');
var turfBBox = require('@turf/bbox').default;
var featureEach = meta.featureEach;
var coordEach = meta.coordEach;
var polygon = helpers.polygon;
var featureCollection = helpers.featureCollection;

/**
 * GeoJSON implementation of [RBush](https://github.com/mourner/rbush#rbush) spatial index.
 *
 * @name rbush
 * @param {number} [maxEntries=9] defines the maximum number of entries in a tree node. 9 (used by default) is a
 * reasonable choice for most applications. Higher value means faster insertion and slower search, and vice versa.
 * @returns {RBush} GeoJSON RBush
 * @example
 * var geojsonRbush = require('geojson-rbush').default;
 * var tree = geojsonRbush();
 */
function geojsonRbush(maxEntries) {
    var tree = rbush(maxEntries);
    /**
     * [insert](https://github.com/mourner/rbush#data-format)
     *
     * @param {Feature} feature insert single GeoJSON Feature
     * @returns {RBush} GeoJSON RBush
     * @example
     * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
     * tree.insert(poly)
     */
    tree.insert = function (feature) {
        if (feature.type !== 'Feature') throw new Error('invalid feature');
        feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
        return rbush.prototype.insert.call(this, feature);
    };

    /**
     * [load](https://github.com/mourner/rbush#bulk-inserting-data)
     *
     * @param {FeatureCollection|Array<Feature>} features load entire GeoJSON FeatureCollection
     * @returns {RBush} GeoJSON RBush
     * @example
     * var polys = turf.polygons([
     *     [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]],
     *     [[[-93, 32], [-83, 32], [-83, 39], [-93, 39], [-93, 32]]]
     * ]);
     * tree.load(polys);
     */
    tree.load = function (features) {
        var load = [];
        // Load an Array of Features
        if (Array.isArray(features)) {
            features.forEach(function (feature) {
                if (feature.type !== 'Feature') throw new Error('invalid features');
                feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
                load.push(feature);
            });
        } else {
            // Load a FeatureCollection
            featureEach(features, function (feature) {
                if (feature.type !== 'Feature') throw new Error('invalid features');
                feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
                load.push(feature);
            });
        }
        return rbush.prototype.load.call(this, load);
    };

    /**
     * [remove](https://github.com/mourner/rbush#removing-data)
     *
     * @param {Feature} feature remove single GeoJSON Feature
     * @param {Function} equals Pass a custom equals function to compare by value for removal.
     * @returns {RBush} GeoJSON RBush
     * @example
     * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
     *
     * tree.remove(poly);
     */
    tree.remove = function (feature, equals) {
        if (feature.type !== 'Feature') throw new Error('invalid feature');
        feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
        return rbush.prototype.remove.call(this, feature, equals);
    };

    /**
     * [clear](https://github.com/mourner/rbush#removing-data)
     *
     * @returns {RBush} GeoJSON Rbush
     * @example
     * tree.clear()
     */
    tree.clear = function () {
        return rbush.prototype.clear.call(this);
    };

    /**
     * [search](https://github.com/mourner/rbush#search)
     *
     * @param {BBox|FeatureCollection|Feature} geojson search with GeoJSON
     * @returns {FeatureCollection} all features that intersects with the given GeoJSON.
     * @example
     * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
     *
     * tree.search(poly);
     */
    tree.search = function (geojson) {
        var features = rbush.prototype.search.call(this, this.toBBox(geojson));
        return featureCollection(features);
    };

    /**
     * [collides](https://github.com/mourner/rbush#collisions)
     *
     * @param {BBox|FeatureCollection|Feature} geojson collides with GeoJSON
     * @returns {boolean} true if there are any items intersecting the given GeoJSON, otherwise false.
     * @example
     * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
     *
     * tree.collides(poly);
     */
    tree.collides = function (geojson) {
        return rbush.prototype.collides.call(this, this.toBBox(geojson));
    };

    /**
     * [all](https://github.com/mourner/rbush#search)
     *
     * @returns {FeatureCollection} all the features in RBush
     * @example
     * tree.all()
     */
    tree.all = function () {
        var features = rbush.prototype.all.call(this);
        return featureCollection(features);
    };

    /**
     * [toJSON](https://github.com/mourner/rbush#export-and-import)
     *
     * @returns {any} export data as JSON object
     * @example
     * var exported = tree.toJSON()
     */
    tree.toJSON = function () {
        return rbush.prototype.toJSON.call(this);
    };

    /**
     * [fromJSON](https://github.com/mourner/rbush#export-and-import)
     *
     * @param {any} json import previously exported data
     * @returns {RBush} GeoJSON RBush
     * @example
     * var exported = {
     *   "children": [
     *     {
     *       "type": "Feature",
     *       "geometry": {
     *         "type": "Point",
     *         "coordinates": [110, 50]
     *       },
     *       "properties": {},
     *       "bbox": [110, 50, 110, 50]
     *     }
     *   ],
     *   "height": 1,
     *   "leaf": true,
     *   "minX": 110,
     *   "minY": 50,
     *   "maxX": 110,
     *   "maxY": 50
     * }
     * tree.fromJSON(exported)
     */
    tree.fromJSON = function (json) {
        return rbush.prototype.fromJSON.call(this, json);
    };

    /**
     * Converts GeoJSON to {minX, minY, maxX, maxY} schema
     *
     * @private
     * @param {BBox|FeatureCollection|Feature} geojson feature(s) to retrieve BBox from
     * @returns {Object} converted to {minX, minY, maxX, maxY}
     */
    tree.toBBox = function (geojson) {
        var bbox;
        if (geojson.bbox) bbox = geojson.bbox;
        else if (Array.isArray(geojson) && geojson.length === 4) bbox = geojson;
        else if (Array.isArray(geojson) && geojson.length === 6) bbox = [geojson[0], geojson[1], geojson[3], geojson[4]];
        else if (geojson.type === 'Feature') bbox = turfBBox(geojson);
        else if (geojson.type === 'FeatureCollection') bbox = turfBBox(geojson);
        else throw new Error('invalid geojson')

        return {
            minX: bbox[0],
            minY: bbox[1],
            maxX: bbox[2],
            maxY: bbox[3]
        };
    };
    return tree;
}

module.exports = geojsonRbush;
module.exports.default = geojsonRbush;

},{"@turf/bbox":3,"@turf/helpers":15,"@turf/meta":52,"rbush":66}],52:[function(require,module,exports){
arguments[4][11][0].apply(exports,arguments)
},{"@turf/helpers":15,"dup":11}],53:[function(require,module,exports){
/**
 * martinez v0.4.3
 * Martinez polygon clipping algorithm, does boolean operation on polygons (multipolygons, polygons with holes etc): intersection, union, difference, xor
 *
 * @author Alex Milevski <info@w8r.name>
 * @license MIT
 * @preserve
 */

(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (factory((global.martinez = {})));
}(this, (function (exports) { 'use strict';

  function DEFAULT_COMPARE (a, b) { return a > b ? 1 : a < b ? -1 : 0; }

  var SplayTree = function SplayTree(compare, noDuplicates) {
    if ( compare === void 0 ) compare = DEFAULT_COMPARE;
    if ( noDuplicates === void 0 ) noDuplicates = false;

    this._compare = compare;
    this._root = null;
    this._size = 0;
    this._noDuplicates = !!noDuplicates;
  };

  var prototypeAccessors = { size: { configurable: true } };


  SplayTree.prototype.rotateLeft = function rotateLeft (x) {
    var y = x.right;
    if (y) {
      x.right = y.left;
      if (y.left) { y.left.parent = x; }
      y.parent = x.parent;
    }

    if (!x.parent)              { this._root = y; }
    else if (x === x.parent.left) { x.parent.left = y; }
    else                        { x.parent.right = y; }
    if (y) { y.left = x; }
    x.parent = y;
  };


  SplayTree.prototype.rotateRight = function rotateRight (x) {
    var y = x.left;
    if (y) {
      x.left = y.right;
      if (y.right) { y.right.parent = x; }
      y.parent = x.parent;
    }

    if (!x.parent)             { this._root = y; }
    else if(x === x.parent.left) { x.parent.left = y; }
    else                       { x.parent.right = y; }
    if (y) { y.right = x; }
    x.parent = y;
  };


  SplayTree.prototype._splay = function _splay (x) {
      var this$1 = this;

    while (x.parent) {
      var p = x.parent;
      if (!p.parent) {
        if (p.left === x) { this$1.rotateRight(p); }
        else            { this$1.rotateLeft(p); }
      } else if (p.left === x && p.parent.left === p) {
        this$1.rotateRight(p.parent);
        this$1.rotateRight(p);
      } else if (p.right === x && p.parent.right === p) {
        this$1.rotateLeft(p.parent);
        this$1.rotateLeft(p);
      } else if (p.left === x && p.parent.right === p) {
        this$1.rotateRight(p);
        this$1.rotateLeft(p);
      } else {
        this$1.rotateLeft(p);
        this$1.rotateRight(p);
      }
    }
  };


  SplayTree.prototype.splay = function splay (x) {
      var this$1 = this;

    var p, gp, ggp, l, r;

    while (x.parent) {
      p = x.parent;
      gp = p.parent;

      if (gp && gp.parent) {
        ggp = gp.parent;
        if (ggp.left === gp) { ggp.left= x; }
        else               { ggp.right = x; }
        x.parent = ggp;
      } else {
        x.parent = null;
        this$1._root = x;
      }

      l = x.left; r = x.right;

      if (x === p.left) { // left
        if (gp) {
          if (gp.left === p) {
            /* zig-zig */
            if (p.right) {
              gp.left = p.right;
              gp.left.parent = gp;
            } else { gp.left = null; }

            p.right = gp;
            gp.parent = p;
          } else {
            /* zig-zag */
            if (l) {
              gp.right = l;
              l.parent = gp;
            } else { gp.right = null; }

            x.left  = gp;
            gp.parent = x;
          }
        }
        if (r) {
          p.left = r;
          r.parent = p;
        } else { p.left = null; }

        x.right= p;
        p.parent = x;
      } else { // right
        if (gp) {
          if (gp.right === p) {
            /* zig-zig */
            if (p.left) {
              gp.right = p.left;
              gp.right.parent = gp;
            } else { gp.right = null; }

            p.left = gp;
            gp.parent = p;
          } else {
            /* zig-zag */
            if (r) {
              gp.left = r;
              r.parent = gp;
            } else { gp.left = null; }

            x.right = gp;
            gp.parent = x;
          }
        }
        if (l) {
          p.right = l;
          l.parent = p;
        } else { p.right = null; }

        x.left = p;
        p.parent = x;
      }
    }
  };


  SplayTree.prototype.replace = function replace (u, v) {
    if (!u.parent) { this._root = v; }
    else if (u === u.parent.left) { u.parent.left = v; }
    else { u.parent.right = v; }
    if (v) { v.parent = u.parent; }
  };


  SplayTree.prototype.minNode = function minNode (u) {
      if ( u === void 0 ) u = this._root;

    if (u) { while (u.left) { u = u.left; } }
    return u;
  };


  SplayTree.prototype.maxNode = function maxNode (u) {
      if ( u === void 0 ) u = this._root;

    if (u) { while (u.right) { u = u.right; } }
    return u;
  };


  SplayTree.prototype.insert = function insert (key, data) {
    var z = this._root;
    var p = null;
    var comp = this._compare;
    var cmp;

    if (this._noDuplicates) {
      while (z) {
        p = z;
        cmp = comp(z.key, key);
        if (cmp === 0) { return; }
        else if (comp(z.key, key) < 0) { z = z.right; }
        else { z = z.left; }
      }
    } else {
      while (z) {
        p = z;
        if (comp(z.key, key) < 0) { z = z.right; }
        else { z = z.left; }
      }
    }

    z = { key: key, data: data, left: null, right: null, parent: p };

    if (!p)                        { this._root = z; }
    else if (comp(p.key, z.key) < 0) { p.right = z; }
    else                           { p.left= z; }

    this.splay(z);
    this._size++;
    return z;
  };


  SplayTree.prototype.find = function find (key) {
    var z  = this._root;
    var comp = this._compare;
    while (z) {
      var cmp = comp(z.key, key);
      if    (cmp < 0) { z = z.right; }
      else if (cmp > 0) { z = z.left; }
      else            { return z; }
    }
    return null;
  };

  /**
   * Whether the tree contains a node with the given key
   * @param{Key} key
   * @return {boolean} true/false
   */
  SplayTree.prototype.contains = function contains (key) {
    var node     = this._root;
    var comparator = this._compare;
    while (node){
      var cmp = comparator(key, node.key);
      if    (cmp === 0) { return true; }
      else if (cmp < 0) { node = node.left; }
      else              { node = node.right; }
    }

    return false;
  };


  SplayTree.prototype.remove = function remove (key) {
    var z = this.find(key);

    if (!z) { return false; }

    this.splay(z);

    if (!z.left) { this.replace(z, z.right); }
    else if (!z.right) { this.replace(z, z.left); }
    else {
      var y = this.minNode(z.right);
      if (y.parent !== z) {
        this.replace(y, y.right);
        y.right = z.right;
        y.right.parent = y;
      }
      this.replace(z, y);
      y.left = z.left;
      y.left.parent = y;
    }

    this._size--;
    return true;
  };


  SplayTree.prototype.removeNode = function removeNode (z) {
    if (!z) { return false; }

    this.splay(z);

    if (!z.left) { this.replace(z, z.right); }
    else if (!z.right) { this.replace(z, z.left); }
    else {
      var y = this.minNode(z.right);
      if (y.parent !== z) {
        this.replace(y, y.right);
        y.right = z.right;
        y.right.parent = y;
      }
      this.replace(z, y);
      y.left = z.left;
      y.left.parent = y;
    }

    this._size--;
    return true;
  };


  SplayTree.prototype.erase = function erase (key) {
    var z = this.find(key);
    if (!z) { return; }

    this.splay(z);

    var s = z.left;
    var t = z.right;

    var sMax = null;
    if (s) {
      s.parent = null;
      sMax = this.maxNode(s);
      this.splay(sMax);
      this._root = sMax;
    }
    if (t) {
      if (s) { sMax.right = t; }
      else { this._root = t; }
      t.parent = sMax;
    }

    this._size--;
  };

  /**
   * Removes and returns the node with smallest key
   * @return {?Node}
   */
  SplayTree.prototype.pop = function pop () {
    var node = this._root, returnValue = null;
    if (node) {
      while (node.left) { node = node.left; }
      returnValue = { key: node.key, data: node.data };
      this.remove(node.key);
    }
    return returnValue;
  };


  /* eslint-disable class-methods-use-this */

  /**
   * Successor node
   * @param{Node} node
   * @return {?Node}
   */
  SplayTree.prototype.next = function next (node) {
    var successor = node;
    if (successor) {
      if (successor.right) {
        successor = successor.right;
        while (successor && successor.left) { successor = successor.left; }
      } else {
        successor = node.parent;
        while (successor && successor.right === node) {
          node = successor; successor = successor.parent;
        }
      }
    }
    return successor;
  };


  /**
   * Predecessor node
   * @param{Node} node
   * @return {?Node}
   */
  SplayTree.prototype.prev = function prev (node) {
    var predecessor = node;
    if (predecessor) {
      if (predecessor.left) {
        predecessor = predecessor.left;
        while (predecessor && predecessor.right) { predecessor = predecessor.right; }
      } else {
        predecessor = node.parent;
        while (predecessor && predecessor.left === node) {
          node = predecessor;
          predecessor = predecessor.parent;
        }
      }
    }
    return predecessor;
  };
  /* eslint-enable class-methods-use-this */


  /**
   * @param{forEachCallback} callback
   * @return {SplayTree}
   */
  SplayTree.prototype.forEach = function forEach (callback) {
    var current = this._root;
    var s = [], done = false, i = 0;

    while (!done) {
      // Reach the left most Node of the current Node
      if (current) {
        // Place pointer to a tree node on the stack
        // before traversing the node's left subtree
        s.push(current);
        current = current.left;
      } else {
        // BackTrack from the empty subtree and visit the Node
        // at the top of the stack; however, if the stack is
        // empty you are done
        if (s.length > 0) {
          current = s.pop();
          callback(current, i++);

          // We have visited the node and its left
          // subtree. Now, it's right subtree's turn
          current = current.right;
        } else { done = true; }
      }
    }
    return this;
  };


  /**
   * Walk key range from `low` to `high`. Stops if `fn` returns a value.
   * @param{Key}    low
   * @param{Key}    high
   * @param{Function} fn
   * @param{*?}     ctx
   * @return {SplayTree}
   */
  SplayTree.prototype.range = function range (low, high, fn, ctx) {
      var this$1 = this;

    var Q = [];
    var compare = this._compare;
    var node = this._root, cmp;

    while (Q.length !== 0 || node) {
      if (node) {
        Q.push(node);
        node = node.left;
      } else {
        node = Q.pop();
        cmp = compare(node.key, high);
        if (cmp > 0) {
          break;
        } else if (compare(node.key, low) >= 0) {
          if (fn.call(ctx, node)) { return this$1; } // stop if smth is returned
        }
        node = node.right;
      }
    }
    return this;
  };

  /**
   * Returns all keys in order
   * @return {Array<Key>}
   */
  SplayTree.prototype.keys = function keys () {
    var current = this._root;
    var s = [], r = [], done = false;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          r.push(current.key);
          current = current.right;
        } else { done = true; }
      }
    }
    return r;
  };


  /**
   * Returns `data` fields of all nodes in order.
   * @return {Array<Value>}
   */
  SplayTree.prototype.values = function values () {
    var current = this._root;
    var s = [], r = [], done = false;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          r.push(current.data);
          current = current.right;
        } else { done = true; }
      }
    }
    return r;
  };


  /**
   * Returns node at given index
   * @param{number} index
   * @return {?Node}
   */
  SplayTree.prototype.at = function at (index) {
    // removed after a consideration, more misleading than useful
    // index = index % this.size;
    // if (index < 0) index = this.size - index;

    var current = this._root;
    var s = [], done = false, i = 0;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          if (i === index) { return current; }
          i++;
          current = current.right;
        } else { done = true; }
      }
    }
    return null;
  };

  /**
   * Bulk-load items. Both array have to be same size
   * @param{Array<Key>}  keys
   * @param{Array<Value>}[values]
   * @param{Boolean}     [presort=false] Pre-sort keys and values, using
   *                                       tree's comparator. Sorting is done
   *                                       in-place
   * @return {AVLTree}
   */
  SplayTree.prototype.load = function load (keys, values, presort) {
      if ( keys === void 0 ) keys = [];
      if ( values === void 0 ) values = [];
      if ( presort === void 0 ) presort = false;

    if (this._size !== 0) { throw new Error('bulk-load: tree is not empty'); }
    var size = keys.length;
    if (presort) { sort(keys, values, 0, size - 1, this._compare); }
    this._root = loadRecursive(null, keys, values, 0, size);
    this._size = size;
    return this;
  };


  SplayTree.prototype.min = function min () {
    var node = this.minNode(this._root);
    if (node) { return node.key; }
    else    { return null; }
  };


  SplayTree.prototype.max = function max () {
    var node = this.maxNode(this._root);
    if (node) { return node.key; }
    else    { return null; }
  };

  SplayTree.prototype.isEmpty = function isEmpty () { return this._root === null; };
  prototypeAccessors.size.get = function () { return this._size; };


  /**
   * Create a tree and load it with items
   * @param{Array<Key>}        keys
   * @param{Array<Value>?}      [values]

   * @param{Function?}          [comparator]
   * @param{Boolean?}           [presort=false] Pre-sort keys and values, using
   *                                             tree's comparator. Sorting is done
   *                                             in-place
   * @param{Boolean?}           [noDuplicates=false] Allow duplicates
   * @return {SplayTree}
   */
  SplayTree.createTree = function createTree (keys, values, comparator, presort, noDuplicates) {
    return new SplayTree(comparator, noDuplicates).load(keys, values, presort);
  };

  Object.defineProperties( SplayTree.prototype, prototypeAccessors );


  function loadRecursive (parent, keys, values, start, end) {
    var size = end - start;
    if (size > 0) {
      var middle = start + Math.floor(size / 2);
      var key    = keys[middle];
      var data   = values[middle];
      var node   = { key: key, data: data, parent: parent };
      node.left    = loadRecursive(node, keys, values, start, middle);
      node.right   = loadRecursive(node, keys, values, middle + 1, end);
      return node;
    }
    return null;
  }


  function sort(keys, values, left, right, compare) {
    if (left >= right) { return; }

    var pivot = keys[(left + right) >> 1];
    var i = left - 1;
    var j = right + 1;

    while (true) {
      do { i++; } while (compare(keys[i], pivot) < 0);
      do { j--; } while (compare(keys[j], pivot) > 0);
      if (i >= j) { break; }

      var tmp = keys[i];
      keys[i] = keys[j];
      keys[j] = tmp;

      tmp = values[i];
      values[i] = values[j];
      values[j] = tmp;
    }

    sort(keys, values,  left,     j, compare);
    sort(keys, values, j + 1, right, compare);
  }

  var NORMAL               = 0;
  var NON_CONTRIBUTING     = 1;
  var SAME_TRANSITION      = 2;
  var DIFFERENT_TRANSITION = 3;

  var INTERSECTION = 0;
  var UNION        = 1;
  var DIFFERENCE   = 2;
  var XOR          = 3;

  /**
   * @param  {SweepEvent} event
   * @param  {SweepEvent} prev
   * @param  {Operation} operation
   */
  function computeFields (event, prev, operation) {
    // compute inOut and otherInOut fields
    if (prev === null) {
      event.inOut      = false;
      event.otherInOut = true;

    // previous line segment in sweepline belongs to the same polygon
    } else {
      if (event.isSubject === prev.isSubject) {
        event.inOut      = !prev.inOut;
        event.otherInOut = prev.otherInOut;

      // previous line segment in sweepline belongs to the clipping polygon
      } else {
        event.inOut      = !prev.otherInOut;
        event.otherInOut = prev.isVertical() ? !prev.inOut : prev.inOut;
      }

      // compute prevInResult field
      if (prev) {
        event.prevInResult = (!inResult(prev, operation) || prev.isVertical())
          ? prev.prevInResult : prev;
      }
    }

    // check if the line segment belongs to the Boolean operation
    event.inResult = inResult(event, operation);
  }


  /* eslint-disable indent */
  function inResult(event, operation) {
    switch (event.type) {
      case NORMAL:
        switch (operation) {
          case INTERSECTION:
            return !event.otherInOut;
          case UNION:
            return event.otherInOut;
          case DIFFERENCE:
            // return (event.isSubject && !event.otherInOut) ||
            //         (!event.isSubject && event.otherInOut);
            return (event.isSubject && event.otherInOut) ||
                    (!event.isSubject && !event.otherInOut);
          case XOR:
            return true;
        }
        break;
      case SAME_TRANSITION:
        return operation === INTERSECTION || operation === UNION;
      case DIFFERENT_TRANSITION:
        return operation === DIFFERENCE;
      case NON_CONTRIBUTING:
        return false;
    }
    return false;
  }
  /* eslint-enable indent */

  var SweepEvent = function SweepEvent (point, left, otherEvent, isSubject, edgeType) {

    /**
     * Is left endpoint?
     * @type {Boolean}
     */
    this.left = left;

    /**
     * @type {Array.<Number>}
     */
    this.point = point;

    /**
     * Other edge reference
     * @type {SweepEvent}
     */
    this.otherEvent = otherEvent;

    /**
     * Belongs to source or clipping polygon
     * @type {Boolean}
     */
    this.isSubject = isSubject;

    /**
     * Edge contribution type
     * @type {Number}
     */
    this.type = edgeType || NORMAL;


    /**
     * In-out transition for the sweepline crossing polygon
     * @type {Boolean}
     */
    this.inOut = false;


    /**
     * @type {Boolean}
     */
    this.otherInOut = false;

    /**
     * Previous event in result?
     * @type {SweepEvent}
     */
    this.prevInResult = null;

    /**
     * Does event belong to result?
     * @type {Boolean}
     */
    this.inResult = false;


    // connection step

    /**
     * @type {Boolean}
     */
    this.resultInOut = false;

    this.isExteriorRing = true;
  };


  /**
   * @param{Array.<Number>}p
   * @return {Boolean}
   */
  SweepEvent.prototype.isBelow = function isBelow (p) {
    var p0 = this.point, p1 = this.otherEvent.point;
    return this.left
      ? (p0[0] - p[0]) * (p1[1] - p[1]) - (p1[0] - p[0]) * (p0[1] - p[1]) > 0
      // signedArea(this.point, this.otherEvent.point, p) > 0 :
      : (p1[0] - p[0]) * (p0[1] - p[1]) - (p0[0] - p[0]) * (p1[1] - p[1]) > 0;
      //signedArea(this.otherEvent.point, this.point, p) > 0;
  };


  /**
   * @param{Array.<Number>}p
   * @return {Boolean}
   */
  SweepEvent.prototype.isAbove = function isAbove (p) {
    return !this.isBelow(p);
  };


  /**
   * @return {Boolean}
   */
  SweepEvent.prototype.isVertical = function isVertical () {
    return this.point[0] === this.otherEvent.point[0];
  };


  SweepEvent.prototype.clone = function clone () {
    var copy = new SweepEvent(
      this.point, this.left, this.otherEvent, this.isSubject, this.type);

    copy.inResult     = this.inResult;
    copy.prevInResult = this.prevInResult;
    copy.isExteriorRing = this.isExteriorRing;
    copy.inOut        = this.inOut;
    copy.otherInOut   = this.otherInOut;

    return copy;
  };

  function equals(p1, p2) {
    if (p1[0] === p2[0]) {
      if (p1[1] === p2[1]) {
        return true;
      } else {
        return false;
      }
    }
    return false;
  }

  // const EPSILON = 1e-9;
  // const abs = Math.abs;
  // TODO https://github.com/w8r/martinez/issues/6#issuecomment-262847164
  // Precision problem.
  //
  // module.exports = function equals(p1, p2) {
  //   return abs(p1[0] - p2[0]) <= EPSILON && abs(p1[1] - p2[1]) <= EPSILON;
  // };

  /**
   * Signed area of the triangle (p0, p1, p2)
   * @param  {Array.<Number>} p0
   * @param  {Array.<Number>} p1
   * @param  {Array.<Number>} p2
   * @return {Number}
   */
  function signedArea(p0, p1, p2) {
    return (p0[0] - p2[0]) * (p1[1] - p2[1]) - (p1[0] - p2[0]) * (p0[1] - p2[1]);
  }

  /**
   * @param  {SweepEvent} e1
   * @param  {SweepEvent} e2
   * @return {Number}
   */
  function compareEvents(e1, e2) {
    var p1 = e1.point;
    var p2 = e2.point;

    // Different x-coordinate
    if (p1[0] > p2[0]) { return 1; }
    if (p1[0] < p2[0]) { return -1; }

    // Different points, but same x-coordinate
    // Event with lower y-coordinate is processed first
    if (p1[1] !== p2[1]) { return p1[1] > p2[1] ? 1 : -1; }

    return specialCases(e1, e2, p1, p2);
  }


  /* eslint-disable no-unused-vars */
  function specialCases(e1, e2, p1, p2) {
    // Same coordinates, but one is a left endpoint and the other is
    // a right endpoint. The right endpoint is processed first
    if (e1.left !== e2.left)
      { return e1.left ? 1 : -1; }

    // const p2 = e1.otherEvent.point, p3 = e2.otherEvent.point;
    // const sa = (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    // Same coordinates, both events
    // are left endpoints or right endpoints.
    // not collinear
    if (signedArea(p1, e1.otherEvent.point, e2.otherEvent.point) !== 0) {
      // the event associate to the bottom segment is processed first
      return (!e1.isBelow(e2.otherEvent.point)) ? 1 : -1;
    }

    return (!e1.isSubject && e2.isSubject) ? 1 : -1;
  }
  /* eslint-enable no-unused-vars */

  /**
   * @param  {SweepEvent} se
   * @param  {Array.<Number>} p
   * @param  {Queue} queue
   * @return {Queue}
   */
  function divideSegment(se, p, queue)  {
    var r = new SweepEvent(p, false, se,            se.isSubject);
    var l = new SweepEvent(p, true,  se.otherEvent, se.isSubject);

    /* eslint-disable no-console */
    if (equals(se.point, se.otherEvent.point)) {

      console.warn('what is that, a collapsed segment?', se);
    }
    /* eslint-enable no-console */

    r.contourId = l.contourId = se.contourId;

    // avoid a rounding error. The left event would be processed after the right event
    if (compareEvents(l, se.otherEvent) > 0) {
      se.otherEvent.left = true;
      l.left = false;
    }

    // avoid a rounding error. The left event would be processed after the right event
    // if (compareEvents(se, r) > 0) {}

    se.otherEvent.otherEvent = l;
    se.otherEvent = r;

    queue.push(l);
    queue.push(r);

    return queue;
  }

  //const EPS = 1e-9;

  /**
   * Finds the magnitude of the cross product of two vectors (if we pretend
   * they're in three dimensions)
   *
   * @param {Object} a First vector
   * @param {Object} b Second vector
   * @private
   * @returns {Number} The magnitude of the cross product
   */
  function crossProduct(a, b) {
    return (a[0] * b[1]) - (a[1] * b[0]);
  }

  /**
   * Finds the dot product of two vectors.
   *
   * @param {Object} a First vector
   * @param {Object} b Second vector
   * @private
   * @returns {Number} The dot product
   */
  function dotProduct(a, b) {
    return (a[0] * b[0]) + (a[1] * b[1]);
  }

  /**
   * Finds the intersection (if any) between two line segments a and b, given the
   * line segments' end points a1, a2 and b1, b2.
   *
   * This algorithm is based on Schneider and Eberly.
   * http://www.cimec.org.ar/~ncalvo/Schneider_Eberly.pdf
   * Page 244.
   *
   * @param {Array.<Number>} a1 point of first line
   * @param {Array.<Number>} a2 point of first line
   * @param {Array.<Number>} b1 point of second line
   * @param {Array.<Number>} b2 point of second line
   * @param {Boolean=}       noEndpointTouch whether to skip single touchpoints
   *                                         (meaning connected segments) as
   *                                         intersections
   * @returns {Array.<Array.<Number>>|Null} If the lines intersect, the point of
   * intersection. If they overlap, the two end points of the overlapping segment.
   * Otherwise, null.
   */
  function intersection (a1, a2, b1, b2, noEndpointTouch) {
    // The algorithm expects our lines in the form P + sd, where P is a point,
    // s is on the interval [0, 1], and d is a vector.
    // We are passed two points. P can be the first point of each pair. The
    // vector, then, could be thought of as the distance (in x and y components)
    // from the first point to the second point.
    // So first, let's make our vectors:
    var va = [a2[0] - a1[0], a2[1] - a1[1]];
    var vb = [b2[0] - b1[0], b2[1] - b1[1]];
    // We also define a function to convert back to regular point form:

    /* eslint-disable arrow-body-style */

    function toPoint(p, s, d) {
      return [
        p[0] + s * d[0],
        p[1] + s * d[1]
      ];
    }

    /* eslint-enable arrow-body-style */

    // The rest is pretty much a straight port of the algorithm.
    var e = [b1[0] - a1[0], b1[1] - a1[1]];
    var kross    = crossProduct(va, vb);
    var sqrKross = kross * kross;
    var sqrLenA  = dotProduct(va, va);
    //const sqrLenB  = dotProduct(vb, vb);

    // Check for line intersection. This works because of the properties of the
    // cross product -- specifically, two vectors are parallel if and only if the
    // cross product is the 0 vector. The full calculation involves relative error
    // to account for possible very small line segments. See Schneider & Eberly
    // for details.
    if (sqrKross > 0/* EPS * sqrLenB * sqLenA */) {
      // If they're not parallel, then (because these are line segments) they
      // still might not actually intersect. This code checks that the
      // intersection point of the lines is actually on both line segments.
      var s = crossProduct(e, vb) / kross;
      if (s < 0 || s > 1) {
        // not on line segment a
        return null;
      }
      var t = crossProduct(e, va) / kross;
      if (t < 0 || t > 1) {
        // not on line segment b
        return null;
      }
      if (s === 0 || s === 1) {
        // on an endpoint of line segment a
        return noEndpointTouch ? null : [toPoint(a1, s, va)];
      }
      if (t === 0 || t === 1) {
        // on an endpoint of line segment b
        return noEndpointTouch ? null : [toPoint(b1, t, vb)];
      }
      return [toPoint(a1, s, va)];
    }

    // If we've reached this point, then the lines are either parallel or the
    // same, but the segments could overlap partially or fully, or not at all.
    // So we need to find the overlap, if any. To do that, we can use e, which is
    // the (vector) difference between the two initial points. If this is parallel
    // with the line itself, then the two lines are the same line, and there will
    // be overlap.
    //const sqrLenE = dotProduct(e, e);
    kross = crossProduct(e, va);
    sqrKross = kross * kross;

    if (sqrKross > 0 /* EPS * sqLenB * sqLenE */) {
    // Lines are just parallel, not the same. No overlap.
      return null;
    }

    var sa = dotProduct(va, e) / sqrLenA;
    var sb = sa + dotProduct(va, vb) / sqrLenA;
    var smin = Math.min(sa, sb);
    var smax = Math.max(sa, sb);

    // this is, essentially, the FindIntersection acting on floats from
    // Schneider & Eberly, just inlined into this function.
    if (smin <= 1 && smax >= 0) {

      // overlap on an end point
      if (smin === 1) {
        return noEndpointTouch ? null : [toPoint(a1, smin > 0 ? smin : 0, va)];
      }

      if (smax === 0) {
        return noEndpointTouch ? null : [toPoint(a1, smax < 1 ? smax : 1, va)];
      }

      if (noEndpointTouch && smin === 0 && smax === 1) { return null; }

      // There's overlap on a segment -- two points of intersection. Return both.
      return [
        toPoint(a1, smin > 0 ? smin : 0, va),
        toPoint(a1, smax < 1 ? smax : 1, va)
      ];
    }

    return null;
  }

  /**
   * @param  {SweepEvent} se1
   * @param  {SweepEvent} se2
   * @param  {Queue}      queue
   * @return {Number}
   */
  function possibleIntersection (se1, se2, queue) {
    // that disallows self-intersecting polygons,
    // did cost us half a day, so I'll leave it
    // out of respect
    // if (se1.isSubject === se2.isSubject) return;
    var inter = intersection(
      se1.point, se1.otherEvent.point,
      se2.point, se2.otherEvent.point
    );

    var nintersections = inter ? inter.length : 0;
    if (nintersections === 0) { return 0; } // no intersection

    // the line segments intersect at an endpoint of both line segments
    if ((nintersections === 1) &&
        (equals(se1.point, se2.point) ||
         equals(se1.otherEvent.point, se2.otherEvent.point))) {
      return 0;
    }

    if (nintersections === 2 && se1.isSubject === se2.isSubject) {
      // if(se1.contourId === se2.contourId){
      // console.warn('Edges of the same polygon overlap',
      //   se1.point, se1.otherEvent.point, se2.point, se2.otherEvent.point);
      // }
      //throw new Error('Edges of the same polygon overlap');
      return 0;
    }

    // The line segments associated to se1 and se2 intersect
    if (nintersections === 1) {

      // if the intersection point is not an endpoint of se1
      if (!equals(se1.point, inter[0]) && !equals(se1.otherEvent.point, inter[0])) {
        divideSegment(se1, inter[0], queue);
      }

      // if the intersection point is not an endpoint of se2
      if (!equals(se2.point, inter[0]) && !equals(se2.otherEvent.point, inter[0])) {
        divideSegment(se2, inter[0], queue);
      }
      return 1;
    }

    // The line segments associated to se1 and se2 overlap
    var events        = [];
    var leftCoincide  = false;
    var rightCoincide = false;

    if (equals(se1.point, se2.point)) {
      leftCoincide = true; // linked
    } else if (compareEvents(se1, se2) === 1) {
      events.push(se2, se1);
    } else {
      events.push(se1, se2);
    }

    if (equals(se1.otherEvent.point, se2.otherEvent.point)) {
      rightCoincide = true;
    } else if (compareEvents(se1.otherEvent, se2.otherEvent) === 1) {
      events.push(se2.otherEvent, se1.otherEvent);
    } else {
      events.push(se1.otherEvent, se2.otherEvent);
    }

    if ((leftCoincide && rightCoincide) || leftCoincide) {
      // both line segments are equal or share the left endpoint
      se2.type = NON_CONTRIBUTING;
      se1.type = (se2.inOut === se1.inOut)
        ? SAME_TRANSITION : DIFFERENT_TRANSITION;

      if (leftCoincide && !rightCoincide) {
        // honestly no idea, but changing events selection from [2, 1]
        // to [0, 1] fixes the overlapping self-intersecting polygons issue
        divideSegment(events[1].otherEvent, events[0].point, queue);
      }
      return 2;
    }

    // the line segments share the right endpoint
    if (rightCoincide) {
      divideSegment(events[0], events[1].point, queue);
      return 3;
    }

    // no line segment includes totally the other one
    if (events[0] !== events[3].otherEvent) {
      divideSegment(events[0], events[1].point, queue);
      divideSegment(events[1], events[2].point, queue);
      return 3;
    }

    // one line segment includes the other one
    divideSegment(events[0], events[1].point, queue);
    divideSegment(events[3].otherEvent, events[2].point, queue);

    return 3;
  }

  /**
   * @param  {SweepEvent} le1
   * @param  {SweepEvent} le2
   * @return {Number}
   */
  function compareSegments(le1, le2) {
    if (le1 === le2) { return 0; }

    // Segments are not collinear
    if (signedArea(le1.point, le1.otherEvent.point, le2.point) !== 0 ||
      signedArea(le1.point, le1.otherEvent.point, le2.otherEvent.point) !== 0) {

      // If they share their left endpoint use the right endpoint to sort
      if (equals(le1.point, le2.point)) { return le1.isBelow(le2.otherEvent.point) ? -1 : 1; }

      // Different left endpoint: use the left endpoint to sort
      if (le1.point[0] === le2.point[0]) { return le1.point[1] < le2.point[1] ? -1 : 1; }

      // has the line segment associated to e1 been inserted
      // into S after the line segment associated to e2 ?
      if (compareEvents(le1, le2) === 1) { return le2.isAbove(le1.point) ? -1 : 1; }

      // The line segment associated to e2 has been inserted
      // into S after the line segment associated to e1
      return le1.isBelow(le2.point) ? -1 : 1;
    }

    if (le1.isSubject === le2.isSubject) { // same polygon
      var p1 = le1.point, p2 = le2.point;
      if (p1[0] === p2[0] && p1[1] === p2[1]/*equals(le1.point, le2.point)*/) {
        p1 = le1.otherEvent.point; p2 = le2.otherEvent.point;
        if (p1[0] === p2[0] && p1[1] === p2[1]) { return 0; }
        else { return le1.contourId > le2.contourId ? 1 : -1; }
      }
    } else { // Segments are collinear, but belong to separate polygons
      return le1.isSubject ? -1 : 1;
    }

    return compareEvents(le1, le2) === 1 ? 1 : -1;
  }

  function subdivide(eventQueue, subject, clipping, sbbox, cbbox, operation) {
    var sweepLine = new SplayTree(compareSegments);
    var sortedEvents = [];

    var rightbound = Math.min(sbbox[2], cbbox[2]);

    var prev, next, begin;

    while (eventQueue.length !== 0) {
      var event = eventQueue.pop();
      sortedEvents.push(event);

      // optimization by bboxes for intersection and difference goes here
      if ((operation === INTERSECTION && event.point[0] > rightbound) ||
          (operation === DIFFERENCE   && event.point[0] > sbbox[2])) {
        break;
      }

      if (event.left) {
        next  = prev = sweepLine.insert(event);
        begin = sweepLine.minNode();

        if (prev !== begin) { prev = sweepLine.prev(prev); }
        else                { prev = null; }

        next = sweepLine.next(next);

        var prevEvent = prev ? prev.key : null;
        var prevprevEvent = (void 0);
        computeFields(event, prevEvent, operation);
        if (next) {
          if (possibleIntersection(event, next.key, eventQueue) === 2) {
            computeFields(event, prevEvent, operation);
            computeFields(event, next.key, operation);
          }
        }

        if (prev) {
          if (possibleIntersection(prev.key, event, eventQueue) === 2) {
            var prevprev = prev;
            if (prevprev !== begin) { prevprev = sweepLine.prev(prevprev); }
            else                    { prevprev = null; }

            prevprevEvent = prevprev ? prevprev.key : null;
            computeFields(prevEvent, prevprevEvent, operation);
            computeFields(event,     prevEvent,     operation);
          }
        }
      } else {
        event = event.otherEvent;
        next = prev = sweepLine.find(event);

        if (prev && next) {

          if (prev !== begin) { prev = sweepLine.prev(prev); }
          else                { prev = null; }

          next = sweepLine.next(next);
          sweepLine.remove(event);

          if (next && prev) {
            possibleIntersection(prev.key, next.key, eventQueue);
          }
        }
      }
    }
    return sortedEvents;
  }

  /**
   * @param  {Array.<SweepEvent>} sortedEvents
   * @return {Array.<SweepEvent>}
   */
  function orderEvents(sortedEvents) {
    var event, i, len, tmp;
    var resultEvents = [];
    for (i = 0, len = sortedEvents.length; i < len; i++) {
      event = sortedEvents[i];
      if ((event.left && event.inResult) ||
        (!event.left && event.otherEvent.inResult)) {
        resultEvents.push(event);
      }
    }
    // Due to overlapping edges the resultEvents array can be not wholly sorted
    var sorted = false;
    while (!sorted) {
      sorted = true;
      for (i = 0, len = resultEvents.length; i < len; i++) {
        if ((i + 1) < len &&
          compareEvents(resultEvents[i], resultEvents[i + 1]) === 1) {
          tmp = resultEvents[i];
          resultEvents[i] = resultEvents[i + 1];
          resultEvents[i + 1] = tmp;
          sorted = false;
        }
      }
    }


    for (i = 0, len = resultEvents.length; i < len; i++) {
      event = resultEvents[i];
      event.pos = i;
    }

    // imagine, the right event is found in the beginning of the queue,
    // when his left counterpart is not marked yet
    for (i = 0, len = resultEvents.length; i < len; i++) {
      event = resultEvents[i];
      if (!event.left) {
        tmp = event.pos;
        event.pos = event.otherEvent.pos;
        event.otherEvent.pos = tmp;
      }
    }

    return resultEvents;
  }


  /**
   * @param  {Number} pos
   * @param  {Array.<SweepEvent>} resultEvents
   * @param  {Object>}    processed
   * @return {Number}
   */
  function nextPos(pos, resultEvents, processed, origIndex) {
    var newPos = pos + 1;
    var length = resultEvents.length;
    if (newPos > length - 1) { return pos - 1; }
    var p  = resultEvents[pos].point;
    var p1 = resultEvents[newPos].point;


    // while in range and not the current one by value
    while (newPos < length && p1[0] === p[0] && p1[1] === p[1]) {
      if (!processed[newPos]) {
        return newPos;
      } else   {
        newPos++;
      }
      p1 = resultEvents[newPos].point;
    }

    newPos = pos - 1;

    while (processed[newPos] && newPos >= origIndex) {
      newPos--;
    }
    return newPos;
  }


  /**
   * @param  {Array.<SweepEvent>} sortedEvents
   * @return {Array.<*>} polygons
   */
  function connectEdges(sortedEvents, operation) {
    var i, len;
    var resultEvents = orderEvents(sortedEvents);

    // "false"-filled array
    var processed = {};
    var result = [];
    var event;

    for (i = 0, len = resultEvents.length; i < len; i++) {
      if (processed[i]) { continue; }
      var contour = [[]];

      if (!resultEvents[i].isExteriorRing) {
        if (operation === DIFFERENCE && !resultEvents[i].isSubject && result.length === 0) {
          result.push(contour);
        } else if (result.length === 0) {
          result.push([[contour]]);
        } else {
          result[result.length - 1].push(contour[0]);
        }
      } else if (operation === DIFFERENCE && !resultEvents[i].isSubject && result.length > 1) {
        result[result.length - 1].push(contour[0]);
      } else {
        result.push(contour);
      }

      var ringId = result.length - 1;
      var pos = i;

      var initial = resultEvents[i].point;
      contour[0].push(initial);

      while (pos >= i) {
        event = resultEvents[pos];
        processed[pos] = true;

        if (event.left) {
          event.resultInOut = false;
          event.contourId   = ringId;
        } else {
          event.otherEvent.resultInOut = true;
          event.otherEvent.contourId   = ringId;
        }

        pos = event.pos;
        processed[pos] = true;
        contour[0].push(resultEvents[pos].point);
        pos = nextPos(pos, resultEvents, processed, i);
      }

      pos = pos === -1 ? i : pos;

      event = resultEvents[pos];
      processed[pos] = processed[event.pos] = true;
      event.otherEvent.resultInOut = true;
      event.otherEvent.contourId   = ringId;
    }

    // Handle if the result is a polygon (eg not multipoly)
    // Commented it again, let's see what do we mean by that
    // if (result.length === 1) result = result[0];
    return result;
  }

  var tinyqueue = TinyQueue;
  var default_1 = TinyQueue;

  function TinyQueue(data, compare) {
      var this$1 = this;

      if (!(this instanceof TinyQueue)) { return new TinyQueue(data, compare); }

      this.data = data || [];
      this.length = this.data.length;
      this.compare = compare || defaultCompare;

      if (this.length > 0) {
          for (var i = (this.length >> 1) - 1; i >= 0; i--) { this$1._down(i); }
      }
  }

  function defaultCompare(a, b) {
      return a < b ? -1 : a > b ? 1 : 0;
  }

  TinyQueue.prototype = {

      push: function (item) {
          this.data.push(item);
          this.length++;
          this._up(this.length - 1);
      },

      pop: function () {
          if (this.length === 0) { return undefined; }

          var top = this.data[0];
          this.length--;

          if (this.length > 0) {
              this.data[0] = this.data[this.length];
              this._down(0);
          }
          this.data.pop();

          return top;
      },

      peek: function () {
          return this.data[0];
      },

      _up: function (pos) {
          var data = this.data;
          var compare = this.compare;
          var item = data[pos];

          while (pos > 0) {
              var parent = (pos - 1) >> 1;
              var current = data[parent];
              if (compare(item, current) >= 0) { break; }
              data[pos] = current;
              pos = parent;
          }

          data[pos] = item;
      },

      _down: function (pos) {
          var this$1 = this;

          var data = this.data;
          var compare = this.compare;
          var halfLength = this.length >> 1;
          var item = data[pos];

          while (pos < halfLength) {
              var left = (pos << 1) + 1;
              var right = left + 1;
              var best = data[left];

              if (right < this$1.length && compare(data[right], best) < 0) {
                  left = right;
                  best = data[right];
              }
              if (compare(best, item) >= 0) { break; }

              data[pos] = best;
              pos = left;
          }

          data[pos] = item;
      }
  };
  tinyqueue.default = default_1;

  var max = Math.max;
  var min = Math.min;

  var contourId = 0;


  function processPolygon(contourOrHole, isSubject, depth, Q, bbox, isExteriorRing) {
    var i, len, s1, s2, e1, e2;
    for (i = 0, len = contourOrHole.length - 1; i < len; i++) {
      s1 = contourOrHole[i];
      s2 = contourOrHole[i + 1];
      e1 = new SweepEvent(s1, false, undefined, isSubject);
      e2 = new SweepEvent(s2, false, e1,        isSubject);
      e1.otherEvent = e2;

      if (s1[0] === s2[0] && s1[1] === s2[1]) {
        continue; // skip collapsed edges, or it breaks
      }

      e1.contourId = e2.contourId = depth;
      if (!isExteriorRing) {
        e1.isExteriorRing = false;
        e2.isExteriorRing = false;
      }
      if (compareEvents(e1, e2) > 0) {
        e2.left = true;
      } else {
        e1.left = true;
      }

      var x = s1[0], y = s1[1];
      bbox[0] = min(bbox[0], x);
      bbox[1] = min(bbox[1], y);
      bbox[2] = max(bbox[2], x);
      bbox[3] = max(bbox[3], y);

      // Pushing it so the queue is sorted from left to right,
      // with object on the left having the highest priority.
      Q.push(e1);
      Q.push(e2);
    }
  }


  function fillQueue(subject, clipping, sbbox, cbbox, operation) {
    var eventQueue = new tinyqueue(null, compareEvents);
    var polygonSet, isExteriorRing, i, ii, j, jj; //, k, kk;

    for (i = 0, ii = subject.length; i < ii; i++) {
      polygonSet = subject[i];
      for (j = 0, jj = polygonSet.length; j < jj; j++) {
        isExteriorRing = j === 0;
        if (isExteriorRing) { contourId++; }
        processPolygon(polygonSet[j], true, contourId, eventQueue, sbbox, isExteriorRing);
      }
    }

    for (i = 0, ii = clipping.length; i < ii; i++) {
      polygonSet = clipping[i];
      for (j = 0, jj = polygonSet.length; j < jj; j++) {
        isExteriorRing = j === 0;
        if (operation === DIFFERENCE) { isExteriorRing = false; }
        if (isExteriorRing) { contourId++; }
        processPolygon(polygonSet[j], false, contourId, eventQueue, cbbox, isExteriorRing);
      }
    }

    return eventQueue;
  }

  var EMPTY = [];


  function trivialOperation(subject, clipping, operation) {
    var result = null;
    if (subject.length * clipping.length === 0) {
      if        (operation === INTERSECTION) {
        result = EMPTY;
      } else if (operation === DIFFERENCE) {
        result = subject;
      } else if (operation === UNION ||
                 operation === XOR) {
        result = (subject.length === 0) ? clipping : subject;
      }
    }
    return result;
  }


  function compareBBoxes(subject, clipping, sbbox, cbbox, operation) {
    var result = null;
    if (sbbox[0] > cbbox[2] ||
        cbbox[0] > sbbox[2] ||
        sbbox[1] > cbbox[3] ||
        cbbox[1] > sbbox[3]) {
      if        (operation === INTERSECTION) {
        result = EMPTY;
      } else if (operation === DIFFERENCE) {
        result = subject;
      } else if (operation === UNION ||
                 operation === XOR) {
        result = subject.concat(clipping);
      }
    }
    return result;
  }


  function boolean(subject, clipping, operation) {
    if (typeof subject[0][0][0] === 'number') {
      subject = [subject];
    }
    if (typeof clipping[0][0][0] === 'number') {
      clipping = [clipping];
    }
    var trivial = trivialOperation(subject, clipping, operation);
    if (trivial) {
      return trivial === EMPTY ? null : trivial;
    }
    var sbbox = [Infinity, Infinity, -Infinity, -Infinity];
    var cbbox = [Infinity, Infinity, -Infinity, -Infinity];

    //console.time('fill queue');
    var eventQueue = fillQueue(subject, clipping, sbbox, cbbox, operation);
    //console.timeEnd('fill queue');

    trivial = compareBBoxes(subject, clipping, sbbox, cbbox, operation);
    if (trivial) {
      return trivial === EMPTY ? null : trivial;
    }
    //console.time('subdivide edges');
    var sortedEvents = subdivide(eventQueue, subject, clipping, sbbox, cbbox, operation);
    //console.timeEnd('subdivide edges');

    //console.time('connect vertices');
    var result = connectEdges(sortedEvents, operation);
    //console.timeEnd('connect vertices');
    return result;
  }

  function union (subject, clipping) {
    return boolean(subject, clipping, UNION);
  }

  function diff (subject, clipping) {
    return boolean(subject, clipping, DIFFERENCE);
  }

  function xor (subject, clipping){
    return boolean(subject, clipping, XOR);
  }

  function intersection$1 (subject, clipping) {
    return boolean(subject, clipping, INTERSECTION);
  }

  /**
   * @enum {Number}
   */
  var operations = { UNION: UNION, DIFFERENCE: DIFFERENCE, INTERSECTION: INTERSECTION, XOR: XOR };

  exports.union = union;
  exports.diff = diff;
  exports.xor = xor;
  exports.intersection = intersection$1;
  exports.operations = operations;

  Object.defineProperty(exports, '__esModule', { value: true });

})));


},{}],54:[function(require,module,exports){
module.exports = function(subject) {
  validateSubject(subject);

  var eventsStorage = createEventsStorage(subject);
  subject.on = eventsStorage.on;
  subject.off = eventsStorage.off;
  subject.fire = eventsStorage.fire;
  return subject;
};

function createEventsStorage(subject) {
  // Store all event listeners to this hash. Key is event name, value is array
  // of callback records.
  //
  // A callback record consists of callback function and its optional context:
  // { 'eventName' => [{callback: function, ctx: object}] }
  var registeredEvents = Object.create(null);

  return {
    on: function (eventName, callback, ctx) {
      if (typeof callback !== 'function') {
        throw new Error('callback is expected to be a function');
      }
      var handlers = registeredEvents[eventName];
      if (!handlers) {
        handlers = registeredEvents[eventName] = [];
      }
      handlers.push({callback: callback, ctx: ctx});

      return subject;
    },

    off: function (eventName, callback) {
      var wantToRemoveAll = (typeof eventName === 'undefined');
      if (wantToRemoveAll) {
        // Killing old events storage should be enough in this case:
        registeredEvents = Object.create(null);
        return subject;
      }

      if (registeredEvents[eventName]) {
        var deleteAllCallbacksForEvent = (typeof callback !== 'function');
        if (deleteAllCallbacksForEvent) {
          delete registeredEvents[eventName];
        } else {
          var callbacks = registeredEvents[eventName];
          for (var i = 0; i < callbacks.length; ++i) {
            if (callbacks[i].callback === callback) {
              callbacks.splice(i, 1);
            }
          }
        }
      }

      return subject;
    },

    fire: function (eventName) {
      var callbacks = registeredEvents[eventName];
      if (!callbacks) {
        return subject;
      }

      var fireArguments;
      if (arguments.length > 1) {
        fireArguments = Array.prototype.splice.call(arguments, 1);
      }
      for(var i = 0; i < callbacks.length; ++i) {
        var callbackInfo = callbacks[i];
        callbackInfo.callback.apply(callbackInfo.ctx, fireArguments);
      }

      return subject;
    }
  };
}

function validateSubject(subject) {
  if (!subject) {
    throw new Error('Eventify cannot use falsy object as events subject');
  }
  var reservedWords = ['on', 'fire', 'off'];
  for (var i = 0; i < reservedWords.length; ++i) {
    if (subject.hasOwnProperty(reservedWords[i])) {
      throw new Error("Subject cannot be eventified, since it already has property '" + reservedWords[i] + "'");
    }
  }
}

},{}],55:[function(require,module,exports){
/**
 * @fileOverview Contains definition of the core graph object.
 */

// TODO: need to change storage layer:
// 1. Be able to get all nodes O(1)
// 2. Be able to get number of links O(1)

/**
 * @example
 *  var graph = require('ngraph.graph')();
 *  graph.addNode(1);     // graph has one node.
 *  graph.addLink(2, 3);  // now graph contains three nodes and one link.
 *
 */
module.exports = createGraph;

var eventify = require('ngraph.events');

/**
 * Creates a new graph
 */
function createGraph(options) {
  // Graph structure is maintained as dictionary of nodes
  // and array of links. Each node has 'links' property which
  // hold all links related to that node. And general links
  // array is used to speed up all links enumeration. This is inefficient
  // in terms of memory, but simplifies coding.
  options = options || {};
  if ('uniqueLinkId' in options) {
    console.warn(
      'ngraph.graph: Starting from version 0.14 `uniqueLinkId` is deprecated.\n' +
      'Use `multigraph` option instead\n',
      '\n',
      'Note: there is also change in default behavior: From now own each graph\n'+
      'is considered to be not a multigraph by default (each edge is unique).'
    );

    options.multigraph = options.uniqueLinkId;
  }

  // Dear reader, the non-multigraphs do not guarantee that there is only
  // one link for a given pair of node. When this option is set to false
  // we can save some memory and CPU (18% faster for non-multigraph);
  if (options.multigraph === undefined) options.multigraph = false;

  var nodes = typeof Object.create === 'function' ? Object.create(null) : {},
    links = [],
    // Hash of multi-edges. Used to track ids of edges between same nodes
    multiEdges = {},
    nodesCount = 0,
    suspendEvents = 0,

    forEachNode = createNodeIterator(),
    createLink = options.multigraph ? createUniqueLink : createSingleLink,

    // Our graph API provides means to listen to graph changes. Users can subscribe
    // to be notified about changes in the graph by using `on` method. However
    // in some cases they don't use it. To avoid unnecessary memory consumption
    // we will not record graph changes until we have at least one subscriber.
    // Code below supports this optimization.
    //
    // Accumulates all changes made during graph updates.
    // Each change element contains:
    //  changeType - one of the strings: 'add', 'remove' or 'update';
    //  node - if change is related to node this property is set to changed graph's node;
    //  link - if change is related to link this property is set to changed graph's link;
    changes = [],
    recordLinkChange = noop,
    recordNodeChange = noop,
    enterModification = noop,
    exitModification = noop;

  // this is our public API:
  var graphPart = {
    /**
     * Adds node to the graph. If node with given id already exists in the graph
     * its data is extended with whatever comes in 'data' argument.
     *
     * @param nodeId the node's identifier. A string or number is preferred.
     * @param [data] additional data for the node being added. If node already
     *   exists its data object is augmented with the new one.
     *
     * @return {node} The newly added node or node with given id if it already exists.
     */
    addNode: addNode,

    /**
     * Adds a link to the graph. The function always create a new
     * link between two nodes. If one of the nodes does not exists
     * a new node is created.
     *
     * @param fromId link start node id;
     * @param toId link end node id;
     * @param [data] additional data to be set on the new link;
     *
     * @return {link} The newly created link
     */
    addLink: addLink,

    /**
     * Removes link from the graph. If link does not exist does nothing.
     *
     * @param link - object returned by addLink() or getLinks() methods.
     *
     * @returns true if link was removed; false otherwise.
     */
    removeLink: removeLink,

    /**
     * Removes node with given id from the graph. If node does not exist in the graph
     * does nothing.
     *
     * @param nodeId node's identifier passed to addNode() function.
     *
     * @returns true if node was removed; false otherwise.
     */
    removeNode: removeNode,

    /**
     * Gets node with given identifier. If node does not exist undefined value is returned.
     *
     * @param nodeId requested node identifier;
     *
     * @return {node} in with requested identifier or undefined if no such node exists.
     */
    getNode: getNode,

    /**
     * Gets number of nodes in this graph.
     *
     * @return number of nodes in the graph.
     */
    getNodesCount: function () {
      return nodesCount;
    },

    /**
     * Gets total number of links in the graph.
     */
    getLinksCount: function () {
      return links.length;
    },

    /**
     * Gets all links (inbound and outbound) from the node with given id.
     * If node with given id is not found null is returned.
     *
     * @param nodeId requested node identifier.
     *
     * @return Array of links from and to requested node if such node exists;
     *   otherwise null is returned.
     */
    getLinks: getLinks,

    /**
     * Invokes callback on each node of the graph.
     *
     * @param {Function(node)} callback Function to be invoked. The function
     *   is passed one argument: visited node.
     */
    forEachNode: forEachNode,

    /**
     * Invokes callback on every linked (adjacent) node to the given one.
     *
     * @param nodeId Identifier of the requested node.
     * @param {Function(node, link)} callback Function to be called on all linked nodes.
     *   The function is passed two parameters: adjacent node and link object itself.
     * @param oriented if true graph treated as oriented.
     */
    forEachLinkedNode: forEachLinkedNode,

    /**
     * Enumerates all links in the graph
     *
     * @param {Function(link)} callback Function to be called on all links in the graph.
     *   The function is passed one parameter: graph's link object.
     *
     * Link object contains at least the following fields:
     *  fromId - node id where link starts;
     *  toId - node id where link ends,
     *  data - additional data passed to graph.addLink() method.
     */
    forEachLink: forEachLink,

    /**
     * Suspend all notifications about graph changes until
     * endUpdate is called.
     */
    beginUpdate: enterModification,

    /**
     * Resumes all notifications about graph changes and fires
     * graph 'changed' event in case there are any pending changes.
     */
    endUpdate: exitModification,

    /**
     * Removes all nodes and links from the graph.
     */
    clear: clear,

    /**
     * Detects whether there is a link between two nodes.
     * Operation complexity is O(n) where n - number of links of a node.
     * NOTE: this function is synonim for getLink()
     *
     * @returns link if there is one. null otherwise.
     */
    hasLink: getLink,

    /**
     * Detects whether there is a node with given id
     * 
     * Operation complexity is O(1)
     * NOTE: this function is synonim for getNode()
     *
     * @returns node if there is one; Falsy value otherwise.
     */
    hasNode: getNode,

    /**
     * Gets an edge between two nodes.
     * Operation complexity is O(n) where n - number of links of a node.
     *
     * @param {string} fromId link start identifier
     * @param {string} toId link end identifier
     *
     * @returns link if there is one. null otherwise.
     */
    getLink: getLink
  };

  // this will add `on()` and `fire()` methods.
  eventify(graphPart);

  monitorSubscribers();

  return graphPart;

  function monitorSubscribers() {
    var realOn = graphPart.on;

    // replace real `on` with our temporary on, which will trigger change
    // modification monitoring:
    graphPart.on = on;

    function on() {
      // now it's time to start tracking stuff:
      graphPart.beginUpdate = enterModification = enterModificationReal;
      graphPart.endUpdate = exitModification = exitModificationReal;
      recordLinkChange = recordLinkChangeReal;
      recordNodeChange = recordNodeChangeReal;

      // this will replace current `on` method with real pub/sub from `eventify`.
      graphPart.on = realOn;
      // delegate to real `on` handler:
      return realOn.apply(graphPart, arguments);
    }
  }

  function recordLinkChangeReal(link, changeType) {
    changes.push({
      link: link,
      changeType: changeType
    });
  }

  function recordNodeChangeReal(node, changeType) {
    changes.push({
      node: node,
      changeType: changeType
    });
  }

  function addNode(nodeId, data) {
    if (nodeId === undefined) {
      throw new Error('Invalid node identifier');
    }

    enterModification();

    var node = getNode(nodeId);
    if (!node) {
      node = new Node(nodeId, data);
      nodesCount++;
      recordNodeChange(node, 'add');
    } else {
      node.data = data;
      recordNodeChange(node, 'update');
    }

    nodes[nodeId] = node;

    exitModification();
    return node;
  }

  function getNode(nodeId) {
    return nodes[nodeId];
  }

  function removeNode(nodeId) {
    var node = getNode(nodeId);
    if (!node) {
      return false;
    }

    enterModification();

    var prevLinks = node.links;
    if (prevLinks) {
      node.links = null;
      for(var i = 0; i < prevLinks.length; ++i) {
        removeLink(prevLinks[i]);
      }
    }

    delete nodes[nodeId];
    nodesCount--;

    recordNodeChange(node, 'remove');

    exitModification();

    return true;
  }


  function addLink(fromId, toId, data) {
    enterModification();

    var fromNode = getNode(fromId) || addNode(fromId);
    var toNode = getNode(toId) || addNode(toId);

    var link = createLink(fromId, toId, data);

    links.push(link);

    // TODO: this is not cool. On large graphs potentially would consume more memory.
    addLinkToNode(fromNode, link);
    if (fromId !== toId) {
      // make sure we are not duplicating links for self-loops
      addLinkToNode(toNode, link);
    }

    recordLinkChange(link, 'add');

    exitModification();

    return link;
  }

  function createSingleLink(fromId, toId, data) {
    var linkId = makeLinkId(fromId, toId);
    return new Link(fromId, toId, data, linkId);
  }

  function createUniqueLink(fromId, toId, data) {
    // TODO: Get rid of this method.
    var linkId = makeLinkId(fromId, toId);
    var isMultiEdge = multiEdges.hasOwnProperty(linkId);
    if (isMultiEdge || getLink(fromId, toId)) {
      if (!isMultiEdge) {
        multiEdges[linkId] = 0;
      }
      var suffix = '@' + (++multiEdges[linkId]);
      linkId = makeLinkId(fromId + suffix, toId + suffix);
    }

    return new Link(fromId, toId, data, linkId);
  }

  function getLinks(nodeId) {
    var node = getNode(nodeId);
    return node ? node.links : null;
  }

  function removeLink(link) {
    if (!link) {
      return false;
    }
    var idx = indexOfElementInArray(link, links);
    if (idx < 0) {
      return false;
    }

    enterModification();

    links.splice(idx, 1);

    var fromNode = getNode(link.fromId);
    var toNode = getNode(link.toId);

    if (fromNode) {
      idx = indexOfElementInArray(link, fromNode.links);
      if (idx >= 0) {
        fromNode.links.splice(idx, 1);
      }
    }

    if (toNode) {
      idx = indexOfElementInArray(link, toNode.links);
      if (idx >= 0) {
        toNode.links.splice(idx, 1);
      }
    }

    recordLinkChange(link, 'remove');

    exitModification();

    return true;
  }

  function getLink(fromNodeId, toNodeId) {
    // TODO: Use sorted links to speed this up
    var node = getNode(fromNodeId),
      i;
    if (!node || !node.links) {
      return null;
    }

    for (i = 0; i < node.links.length; ++i) {
      var link = node.links[i];
      if (link.fromId === fromNodeId && link.toId === toNodeId) {
        return link;
      }
    }

    return null; // no link.
  }

  function clear() {
    enterModification();
    forEachNode(function(node) {
      removeNode(node.id);
    });
    exitModification();
  }

  function forEachLink(callback) {
    var i, length;
    if (typeof callback === 'function') {
      for (i = 0, length = links.length; i < length; ++i) {
        callback(links[i]);
      }
    }
  }

  function forEachLinkedNode(nodeId, callback, oriented) {
    var node = getNode(nodeId);

    if (node && node.links && typeof callback === 'function') {
      if (oriented) {
        return forEachOrientedLink(node.links, nodeId, callback);
      } else {
        return forEachNonOrientedLink(node.links, nodeId, callback);
      }
    }
  }

  function forEachNonOrientedLink(links, nodeId, callback) {
    var quitFast;
    for (var i = 0; i < links.length; ++i) {
      var link = links[i];
      var linkedNodeId = link.fromId === nodeId ? link.toId : link.fromId;

      quitFast = callback(nodes[linkedNodeId], link);
      if (quitFast) {
        return true; // Client does not need more iterations. Break now.
      }
    }
  }

  function forEachOrientedLink(links, nodeId, callback) {
    var quitFast;
    for (var i = 0; i < links.length; ++i) {
      var link = links[i];
      if (link.fromId === nodeId) {
        quitFast = callback(nodes[link.toId], link);
        if (quitFast) {
          return true; // Client does not need more iterations. Break now.
        }
      }
    }
  }

  // we will not fire anything until users of this library explicitly call `on()`
  // method.
  function noop() {}

  // Enter, Exit modification allows bulk graph updates without firing events.
  function enterModificationReal() {
    suspendEvents += 1;
  }

  function exitModificationReal() {
    suspendEvents -= 1;
    if (suspendEvents === 0 && changes.length > 0) {
      graphPart.fire('changed', changes);
      changes.length = 0;
    }
  }

  function createNodeIterator() {
    // Object.keys iterator is 1.3x faster than `for in` loop.
    // See `https://github.com/anvaka/ngraph.graph/tree/bench-for-in-vs-obj-keys`
    // branch for perf test
    return Object.keys ? objectKeysIterator : forInIterator;
  }

  function objectKeysIterator(callback) {
    if (typeof callback !== 'function') {
      return;
    }

    var keys = Object.keys(nodes);
    for (var i = 0; i < keys.length; ++i) {
      if (callback(nodes[keys[i]])) {
        return true; // client doesn't want to proceed. Return.
      }
    }
  }

  function forInIterator(callback) {
    if (typeof callback !== 'function') {
      return;
    }
    var node;

    for (node in nodes) {
      if (callback(nodes[node])) {
        return true; // client doesn't want to proceed. Return.
      }
    }
  }
}

// need this for old browsers. Should this be a separate module?
function indexOfElementInArray(element, array) {
  if (!array) return -1;

  if (array.indexOf) {
    return array.indexOf(element);
  }

  var len = array.length,
    i;

  for (i = 0; i < len; i += 1) {
    if (array[i] === element) {
      return i;
    }
  }

  return -1;
}

/**
 * Internal structure to represent node;
 */
function Node(id, data) {
  this.id = id;
  this.links = null;
  this.data = data;
}

function addLinkToNode(node, link) {
  if (node.links) {
    node.links.push(link);
  } else {
    node.links = [link];
  }
}

/**
 * Internal structure to represent links;
 */
function Link(fromId, toId, data, id) {
  this.fromId = fromId;
  this.toId = toId;
  this.data = data;
  this.id = id;
}

function hashCode(str) {
  var hash = 0, i, chr, len;
  if (str.length == 0) return hash;
  for (i = 0, len = str.length; i < len; i++) {
    chr   = str.charCodeAt(i);
    hash  = ((hash << 5) - hash) + chr;
    hash |= 0; // Convert to 32bit integer
  }
  return hash;
}

function makeLinkId(fromId, toId) {
  return fromId.toString() + '👉 ' + toId.toString();
}

},{"ngraph.events":54}],56:[function(require,module,exports){
/**
 * Based on https://github.com/mourner/tinyqueue
 * Copyright (c) 2017, Vladimir Agafonkin https://github.com/mourner/tinyqueue/blob/master/LICENSE
 * 
 * Adapted for PathFinding needs by @anvaka
 * Copyright (c) 2017, Andrei Kashcha
 */
module.exports = NodeHeap;

function NodeHeap(data, options) {
  if (!(this instanceof NodeHeap)) return new NodeHeap(data, options);

  if (!Array.isArray(data)) {
    // assume first argument is our config object;
    options = data;
    data = [];
  }

  options = options || {};

  this.data = data || [];
  this.length = this.data.length;
  this.compare = options.compare || defaultCompare;
  this.setNodeId = options.setNodeId || noop;

  if (this.length > 0) {
    for (var i = (this.length >> 1); i >= 0; i--) this._down(i);
  }

  if (options.setNodeId) {
    for (var i = 0; i < this.length; ++i) {
      this.setNodeId(this.data[i], i);
    }
  }
}

function noop() {}

function defaultCompare(a, b) {
  return a - b;
}

NodeHeap.prototype = {

  push: function (item) {
    this.data.push(item);
    this.setNodeId(item, this.length);
    this.length++;
    this._up(this.length - 1);
  },

  pop: function () {
    if (this.length === 0) return undefined;

    var top = this.data[0];
    this.length--;

    if (this.length > 0) {
      this.data[0] = this.data[this.length];
      this.setNodeId(this.data[0], 0);
      this._down(0);
    }
    this.data.pop();

    return top;
  },

  peek: function () {
    return this.data[0];
  },

  updateItem: function (pos) {
    this._down(pos);
    this._up(pos);
  },

  _up: function (pos) {
    var data = this.data;
    var compare = this.compare;
    var setNodeId = this.setNodeId;
    var item = data[pos];

    while (pos > 0) {
      var parent = (pos - 1) >> 1;
      var current = data[parent];
      if (compare(item, current) >= 0) break;
        data[pos] = current;

       setNodeId(current, pos);
       pos = parent;
    }

    data[pos] = item;
    setNodeId(item, pos);
  },

  _down: function (pos) {
    var data = this.data;
    var compare = this.compare;
    var halfLength = this.length >> 1;
    var item = data[pos];
    var setNodeId = this.setNodeId;

    while (pos < halfLength) {
      var left = (pos << 1) + 1;
      var right = left + 1;
      var best = data[left];

      if (right < this.length && compare(data[right], best) < 0) {
        left = right;
        best = data[right];
      }
      if (compare(best, item) >= 0) break;

      data[pos] = best;
      setNodeId(best, pos);
      pos = left;
    }

    data[pos] = item;
    setNodeId(item, pos);
  }
};
},{}],57:[function(require,module,exports){
/**
 * Performs suboptimal, greed A Star path finding.
 * This finder does not necessary finds the shortest path. The path
 * that it finds is very close to the shortest one. It is very fast though.
 */
module.exports = aStarBi;

var NodeHeap = require('./NodeHeap');
var makeSearchStatePool = require('./makeSearchStatePool');
var heuristics = require('./heuristics');
var defaultSettings = require('./defaultSettings');

var BY_FROM = 1;
var BY_TO = 2;
var NO_PATH = defaultSettings.NO_PATH;

module.exports.l2 = heuristics.l2;
module.exports.l1 = heuristics.l1;

/**
 * Creates a new instance of pathfinder. A pathfinder has just one method:
 * `find(fromId, toId)`, it may be extended in future.
 * 
 * NOTE: Algorithm implemented in this code DOES NOT find optimal path.
 * Yet the path that it finds is always near optimal, and it finds it very fast.
 * 
 * @param {ngraph.graph} graph instance. See https://github.com/anvaka/ngraph.graph
 * 
 * @param {Object} options that configures search
 * @param {Function(a, b)} options.heuristic - a function that returns estimated distance between
 * nodes `a` and `b`.  Defaults function returns 0, which makes this search equivalent to Dijkstra search.
 * @param {Function(a, b)} options.distance - a function that returns actual distance between two
 * nodes `a` and `b`. By default this is set to return graph-theoretical distance (always 1);
 * 
 * @returns {Object} A pathfinder with single method `find()`.
 */
function aStarBi(graph, options) {
  options = options || {};
  // whether traversal should be considered over oriented graph.
  var oriented = options.oriented;

  var heuristic = options.heuristic;
  if (!heuristic) heuristic = defaultSettings.heuristic;

  var distance = options.distance;
  if (!distance) distance = defaultSettings.distance;
  var pool = makeSearchStatePool();

  return {
    find: find
  };

  function find(fromId, toId) {
    // Not sure if we should return NO_PATH or throw. Throw seem to be more
    // helpful to debug errors. So, throwing.
    var from = graph.getNode(fromId);
    if (!from) throw new Error('fromId is not defined in this graph: ' + fromId);
    var to = graph.getNode(toId);
    if (!to) throw new Error('toId is not defined in this graph: ' + toId);

    if (from === to) return [from]; // trivial case.

    pool.reset();

    var callVisitor = oriented ? orientedVisitor : nonOrientedVisitor;

    // Maps nodeId to NodeSearchState.
    var nodeState = new Map();

    var openSetFrom = new NodeHeap({
      compare: defaultSettings.compareFScore,
      setNodeId: defaultSettings.setHeapIndex
    });

    var openSetTo = new NodeHeap({
      compare: defaultSettings.compareFScore,
      setNodeId: defaultSettings.setHeapIndex
    });


    var startNode = pool.createNewState(from);
    nodeState.set(fromId, startNode);

    // For the first node, fScore is completely heuristic.
    startNode.fScore = heuristic(from, to);
    // The cost of going from start to start is zero.
    startNode.distanceToSource = 0;
    openSetFrom.push(startNode);
    startNode.open = BY_FROM;

    var endNode = pool.createNewState(to);
    endNode.fScore = heuristic(to, from);
    endNode.distanceToSource = 0;
    openSetTo.push(endNode);
    endNode.open = BY_TO;

    // Cost of the best solution found so far. Used for accurate termination
    var lMin = Number.POSITIVE_INFINITY;
    var minFrom;
    var minTo;

    var currentSet = openSetFrom;
    var currentOpener = BY_FROM;

    while (openSetFrom.length > 0 && openSetTo.length > 0) {
      if (openSetFrom.length < openSetTo.length) {
        // we pick a set with less elements
        currentOpener = BY_FROM;
        currentSet = openSetFrom;
      } else {
        currentOpener = BY_TO;
        currentSet = openSetTo;
      }

      var current = currentSet.pop();

      // no need to visit this node anymore
      current.closed = true;

      if (current.distanceToSource > lMin) continue;

      graph.forEachLinkedNode(current.node.id, callVisitor);

      if (minFrom && minTo) {
        // This is not necessary the best path, but we are so greedy that we
        // can't resist:
        return reconstructBiDirectionalPath(minFrom, minTo);
      }
    }

    return NO_PATH; // No path.

    function nonOrientedVisitor(otherNode, link) {
      return visitNode(otherNode, link, current);
    }

    function orientedVisitor(otherNode, link) {
      // For oritned graphs we need to reverse graph, when traveling
      // backwards. So, we use non-oriented ngraph's traversal, and 
      // filter link orientation here.
      if (currentOpener === BY_FROM) {
        if (link.fromId === current.node.id) return visitNode(otherNode, link, current)
      } else if (currentOpener === BY_TO) {
        if (link.toId === current.node.id) return visitNode(otherNode, link, current);
      }
    }

    function canExit(currentNode) {
      var opener = currentNode.open
      if (opener && opener !== currentOpener) {
        return true;
      }

      return false;
    }

    function reconstructBiDirectionalPath(a, b) {
      var pathOfNodes = [];
      var aParent = a;
      while(aParent) {
        pathOfNodes.push(aParent.node);
        aParent = aParent.parent;
      }
      var bParent = b;
      while (bParent) {
        pathOfNodes.unshift(bParent.node);
        bParent = bParent.parent
      }
      return pathOfNodes;
    }

    function visitNode(otherNode, link, cameFrom) {
      var otherSearchState = nodeState.get(otherNode.id);
      if (!otherSearchState) {
        otherSearchState = pool.createNewState(otherNode);
        nodeState.set(otherNode.id, otherSearchState);
      }

      if (otherSearchState.closed) {
        // Already processed this node.
        return;
      }

      if (canExit(otherSearchState, cameFrom)) {
        // this node was opened by alternative opener. The sets intersect now,
        // we found an optimal path, that goes through *this* node. However, there
        // is no guarantee that this is the global optimal solution path.

        var potentialLMin = otherSearchState.distanceToSource + cameFrom.distanceToSource;
        if (potentialLMin < lMin) {
          minFrom = otherSearchState;
          minTo = cameFrom
          lMin = potentialLMin;
        }
        // we are done with this node.
        return;
      }

      var tentativeDistance = cameFrom.distanceToSource + distance(otherSearchState.node, cameFrom.node, link);

      if (tentativeDistance >= otherSearchState.distanceToSource) {
        // This would only make our path longer. Ignore this route.
        return;
      }

      // Choose target based on current working set:
      var target = (currentOpener === BY_FROM) ? to : from;
      var newFScore = tentativeDistance + heuristic(otherSearchState.node, target);
      if (newFScore >= lMin) {
        // this can't be optimal path, as we have already found a shorter path.
        return;
      }
      otherSearchState.fScore = newFScore;

      if (otherSearchState.open === 0) {
        // Remember this node in the current set
        currentSet.push(otherSearchState);
        currentSet.updateItem(otherSearchState.heapIndex);

        otherSearchState.open = currentOpener;
      }

      // bingo! we found shorter path:
      otherSearchState.parent = cameFrom;
      otherSearchState.distanceToSource = tentativeDistance;
    }
  }
}

},{"./NodeHeap":56,"./defaultSettings":59,"./heuristics":60,"./makeSearchStatePool":61}],58:[function(require,module,exports){
/**
 * Performs a uni-directional A Star search on graph.
 * 
 * We will try to minimize f(n) = g(n) + h(n), where
 * g(n) is actual distance from source node to `n`, and
 * h(n) is heuristic distance from `n` to target node.
 */
module.exports = aStarPathSearch;

var NodeHeap = require('./NodeHeap');
var makeSearchStatePool = require('./makeSearchStatePool');
var heuristics = require('./heuristics');
var defaultSettings = require('./defaultSettings.js');

var NO_PATH = defaultSettings.NO_PATH;

module.exports.l2 = heuristics.l2;
module.exports.l1 = heuristics.l1;

/**
 * Creates a new instance of pathfinder. A pathfinder has just one method:
 * `find(fromId, toId)`, it may be extended in future.
 * 
 * @param {ngraph.graph} graph instance. See https://github.com/anvaka/ngraph.graph
 * @param {Object} options that configures search
 * @param {Function(a, b)} options.heuristic - a function that returns estimated distance between
 * nodes `a` and `b`. This function should never overestimate actual distance between two
 * nodes (otherwise the found path will not be the shortest). Defaults function returns 0,
 * which makes this search equivalent to Dijkstra search.
 * @param {Function(a, b)} options.distance - a function that returns actual distance between two
 * nodes `a` and `b`. By default this is set to return graph-theoretical distance (always 1);
 * 
 * @returns {Object} A pathfinder with single method `find()`.
 */
function aStarPathSearch(graph, options) {
  options = options || {};
  // whether traversal should be considered over oriented graph.
  var oriented = options.oriented;

  var heuristic = options.heuristic;
  if (!heuristic) heuristic = defaultSettings.heuristic;

  var distance = options.distance;
  if (!distance) distance = defaultSettings.distance;
  var pool = makeSearchStatePool();

  return {
    /**
     * Finds a path between node `fromId` and `toId`.
     * @returns {Array} of nodes between `toId` and `fromId`. Empty array is returned
     * if no path is found.
     */
    find: find
  };

  function find(fromId, toId) {
    var from = graph.getNode(fromId);
    if (!from) throw new Error('fromId is not defined in this graph: ' + fromId);
    var to = graph.getNode(toId);
    if (!to) throw new Error('toId is not defined in this graph: ' + toId);
    pool.reset();

    // Maps nodeId to NodeSearchState.
    var nodeState = new Map();

    // the nodes that we still need to evaluate
    var openSet = new NodeHeap({
      compare: defaultSettings.compareFScore,
      setNodeId: defaultSettings.setHeapIndex
    });

    var startNode = pool.createNewState(from);
    nodeState.set(fromId, startNode);

    // For the first node, fScore is completely heuristic.
    startNode.fScore = heuristic(from, to);

    // The cost of going from start to start is zero.
    startNode.distanceToSource = 0;
    openSet.push(startNode);
    startNode.open = 1;

    var cameFrom;

    while (openSet.length > 0) {
      cameFrom = openSet.pop();
      if (goalReached(cameFrom, to)) return reconstructPath(cameFrom);

      // no need to visit this node anymore
      cameFrom.closed = true;
      graph.forEachLinkedNode(cameFrom.node.id, visitNeighbour, oriented);
    }

    // If we got here, then there is no path.
    return NO_PATH;

    function visitNeighbour(otherNode, link) {
      var otherSearchState = nodeState.get(otherNode.id);
      if (!otherSearchState) {
        otherSearchState = pool.createNewState(otherNode);
        nodeState.set(otherNode.id, otherSearchState);
      }

      if (otherSearchState.closed) {
        // Already processed this node.
        return;
      }
      if (otherSearchState.open === 0) {
        // Remember this node.
        openSet.push(otherSearchState);
        otherSearchState.open = 1;
      }

      var tentativeDistance = cameFrom.distanceToSource + distance(otherNode, cameFrom.node, link);
      if (tentativeDistance >= otherSearchState.distanceToSource) {
        // This would only make our path longer. Ignore this route.
        return;
      }

      // bingo! we found shorter path:
      otherSearchState.parent = cameFrom;
      otherSearchState.distanceToSource = tentativeDistance;
      otherSearchState.fScore = tentativeDistance + heuristic(otherSearchState.node, to);

      openSet.updateItem(otherSearchState.heapIndex);
    }
  }
}

function goalReached(searchState, targetNode) {
  return searchState.node === targetNode;
}

function reconstructPath(searchState) {
  var path = [searchState.node];
  var parent = searchState.parent;

  while (parent) {
    path.push(parent.node);
    parent = parent.parent;
  }

  return path;
}

},{"./NodeHeap":56,"./defaultSettings.js":59,"./heuristics":60,"./makeSearchStatePool":61}],59:[function(require,module,exports){
// We reuse instance of array, but we trie to freeze it as well,
// so that consumers don't modify it. Maybe it's a bad idea.
var NO_PATH = [];
if (typeof Object.freeze === 'function') Object.freeze(NO_PATH);

module.exports = {
  // Path search settings
  heuristic: blindHeuristic,
  distance: constantDistance,
  compareFScore: compareFScore,
  NO_PATH: NO_PATH,

  // heap settings
  setHeapIndex: setHeapIndex,

  // nba:
  setH1: setH1,
  setH2: setH2,
  compareF1Score: compareF1Score,
  compareF2Score: compareF2Score,
}

function blindHeuristic(/* a, b */) {
  // blind heuristic makes this search equal to plain Dijkstra path search.
  return 0;
}

function constantDistance(/* a, b */) {
  return 1;
}

function compareFScore(a, b) {
  var result = a.fScore - b.fScore;
  // TODO: Can I improve speed with smarter ties-breaking?
  // I tried distanceToSource, but it didn't seem to have much effect
  return result;
}

function setHeapIndex(nodeSearchState, heapIndex) {
  nodeSearchState.heapIndex = heapIndex;
}

function compareF1Score(a, b) {
  return a.f1 - b.f1;
}

function compareF2Score(a, b) {
  return a.f2 - b.f2;
}

function setH1(node, heapIndex) {
  node.h1 = heapIndex;
}

function setH2(node, heapIndex) {
  node.h2 = heapIndex;
}
},{}],60:[function(require,module,exports){
module.exports = {
  l2: l2,
  l1: l1
};

/**
 * Euclid distance (l2 norm);
 * 
 * @param {*} a 
 * @param {*} b 
 */
function l2(a, b) {
  var dx = a.x - b.x;
  var dy = a.y - b.y;
  return Math.sqrt(dx * dx + dy * dy);
}

/**
 * Manhattan distance (l1 norm);
 * @param {*} a 
 * @param {*} b 
 */
function l1(a, b) {
  var dx = a.x - b.x;
  var dy = a.y - b.y;
  return Math.abs(dx) + Math.abs(dy);
}

},{}],61:[function(require,module,exports){
/**
 * This class represents a single search node in the exploration tree for
 * A* algorithm.
 * 
 * @param {Object} node  original node in the graph
 */
function NodeSearchState(node) {
  this.node = node;

  // How we came to this node?
  this.parent = null;

  this.closed = false;
  this.open = 0;

  this.distanceToSource = Number.POSITIVE_INFINITY;
  // the f(n) = g(n) + h(n) value
  this.fScore = Number.POSITIVE_INFINITY;

  // used to reconstruct heap when fScore is updated.
  this.heapIndex = -1;
};

function makeSearchStatePool() {
  var currentInCache = 0;
  var nodeCache = [];

  return {
    createNewState: createNewState,
    reset: reset
  };

  function reset() {
    currentInCache = 0;
  }

  function createNewState(node) {
    var cached = nodeCache[currentInCache];
    if (cached) {
      // TODO: This almost duplicates constructor code. Not sure if
      // it would impact performance if I move this code into a function
      cached.node = node;
      // How we came to this node?
      cached.parent = null;

      cached.closed = false;
      cached.open = 0;

      cached.distanceToSource = Number.POSITIVE_INFINITY;
      // the f(n) = g(n) + h(n) value
      cached.fScore = Number.POSITIVE_INFINITY;

      // used to reconstruct heap when fScore is updated.
      cached.heapIndex = -1;

    } else {
      cached = new NodeSearchState(node);
      nodeCache[currentInCache] = cached;
    }
    currentInCache++;
    return cached;
  }
}
module.exports = makeSearchStatePool;
},{}],62:[function(require,module,exports){
module.exports = nba;

var NodeHeap = require('../NodeHeap');
var heuristics = require('../heuristics');
var defaultSettings = require('../defaultSettings.js');
var makeNBASearchStatePool = require('./makeNBASearchStatePool.js');

var NO_PATH = defaultSettings.NO_PATH;

module.exports.l2 = heuristics.l2;
module.exports.l1 = heuristics.l1;

/**
 * Creates a new instance of pathfinder. A pathfinder has just one method:
 * `find(fromId, toId)`.
 * 
 * This is implementation of the NBA* algorithm described in 
 * 
 *  "Yet another bidirectional algorithm for shortest paths" paper by Wim Pijls and Henk Post
 * 
 * The paper is available here: https://repub.eur.nl/pub/16100/ei2009-10.pdf
 * 
 * @param {ngraph.graph} graph instance. See https://github.com/anvaka/ngraph.graph
 * @param {Object} options that configures search
 * @param {Function(a, b)} options.heuristic - a function that returns estimated distance between
 * nodes `a` and `b`. This function should never overestimate actual distance between two
 * nodes (otherwise the found path will not be the shortest). Defaults function returns 0,
 * which makes this search equivalent to Dijkstra search.
 * @param {Function(a, b)} options.distance - a function that returns actual distance between two
 * nodes `a` and `b`. By default this is set to return graph-theoretical distance (always 1);
 * 
 * @returns {Object} A pathfinder with single method `find()`.
 */
function nba(graph, options) {
  options = options || {};
  // whether traversal should be considered over oriented graph.
  var oriented = options.oriented;
  var quitFast = options.quitFast;

  var heuristic = options.heuristic;
  if (!heuristic) heuristic = defaultSettings.heuristic;

  var distance = options.distance;
  if (!distance) distance = defaultSettings.distance;

  // During stress tests I noticed that garbage collection was one of the heaviest
  // contributors to the algorithm's speed. So I'm using an object pool to recycle nodes.
  var pool = makeNBASearchStatePool();

  return {
    /**
     * Finds a path between node `fromId` and `toId`.
     * @returns {Array} of nodes between `toId` and `fromId`. Empty array is returned
     * if no path is found.
     */
    find: find
  };

  function find(fromId, toId) {
    // I must apologize for the code duplication. This was the easiest way for me to
    // implement the algorithm fast.
    var from = graph.getNode(fromId);
    if (!from) throw new Error('fromId is not defined in this graph: ' + fromId);
    var to = graph.getNode(toId);
    if (!to) throw new Error('toId is not defined in this graph: ' + toId);

    pool.reset();

    // I must also apologize for somewhat cryptic names. The NBA* is bi-directional
    // search algorithm, which means it runs two searches in parallel. One runs
    // from source node to target, while the other one runs from target to source.
    // Everywhere where you see `1` it means it's for the forward search. `2` is for 
    // backward search.

    // For oriented graph path finding, we need to reverse the graph, so that
    // backward search visits correct link. Obviously we don't want to duplicate
    // the graph, instead we always traverse the graph as non-oriented, and filter
    // edges in `visitN1Oriented/visitN2Oritented`
    var forwardVisitor = oriented ? visitN1Oriented : visitN1;
    var reverseVisitor = oriented ? visitN2Oriented : visitN2;

    // Maps nodeId to NBASearchState.
    var nodeState = new Map();

    // These two heaps store nodes by their underestimated values.
    var open1Set = new NodeHeap({
      compare: defaultSettings.compareF1Score,
      setNodeId: defaultSettings.setH1
    });
    var open2Set = new NodeHeap({
      compare: defaultSettings.compareF2Score,
      setNodeId: defaultSettings.setH2
    });

    // This is where both searches will meet.
    var minNode;

    // The smallest path length seen so far is stored here:
    var lMin = Number.POSITIVE_INFINITY;

    // We start by putting start/end nodes to the corresponding heaps
    var startNode = pool.createNewState(from);
    nodeState.set(fromId, startNode); 
    startNode.g1 = 0;
    var f1 = heuristic(from, to);
    startNode.f1 = f1;
    open1Set.push(startNode);

    var endNode = pool.createNewState(to);
    nodeState.set(toId, endNode);
    endNode.g2 = 0;
    var f2 = f1; // they should agree originally
    endNode.f2 = f2;
    open2Set.push(endNode)

    // the `cameFrom` variable is accessed by both searches, so that we can store parents.
    var cameFrom;

    // this is the main algorithm loop:
    while (open2Set.length && open1Set.length) {
      if (open1Set.length < open2Set.length) {
        forwardSearch();
      } else {
        reverseSearch();
      }

      if (quitFast && minNode) break;
    }

    // If we got here, then there is no path.
    var path = reconstructPath(minNode);
    return path; // the public API is over

    function forwardSearch() {
      cameFrom = open1Set.pop();
      if (cameFrom.closed) {
        return;
      }

      cameFrom.closed = true;

      if (cameFrom.f1 < lMin && (cameFrom.g1 + f2 - heuristic(from, cameFrom.node)) < lMin) {
        graph.forEachLinkedNode(cameFrom.node.id, forwardVisitor);
      }

      if (open1Set.length > 0) {
        f1 = open1Set.peek().f1;
      } 
    }

    function reverseSearch() {
      cameFrom = open2Set.pop();
      if (cameFrom.closed) {
        return;
      }
      cameFrom.closed = true;

      if (cameFrom.f2 < lMin && (cameFrom.g2 + f1 - heuristic(cameFrom.node, to)) < lMin) {
        graph.forEachLinkedNode(cameFrom.node.id, reverseVisitor);
      }

      if (open2Set.length > 0) {
        f2 = open2Set.peek().f2;
      }
    }

    function visitN1(otherNode, link) {
      var otherSearchState = nodeState.get(otherNode.id);
      if (!otherSearchState) {
        otherSearchState = pool.createNewState(otherNode);
        nodeState.set(otherNode.id, otherSearchState);
      }

      if (otherSearchState.closed) return;

      var tentativeDistance = cameFrom.g1 + distance(cameFrom.node, otherNode, link);

      if (tentativeDistance < otherSearchState.g1) {
        otherSearchState.g1 = tentativeDistance;
        otherSearchState.f1 = tentativeDistance + heuristic(otherSearchState.node, to);
        otherSearchState.p1 = cameFrom;
        if (otherSearchState.h1 < 0) {
          open1Set.push(otherSearchState);
        } else {
          open1Set.updateItem(otherSearchState.h1);
        }
      }
      var potentialMin = otherSearchState.g1 + otherSearchState.g2;
      if (potentialMin < lMin) { 
        lMin = potentialMin;
        minNode = otherSearchState;
      }
    }

    function visitN2(otherNode, link) {
      var otherSearchState = nodeState.get(otherNode.id);
      if (!otherSearchState) {
        otherSearchState = pool.createNewState(otherNode);
        nodeState.set(otherNode.id, otherSearchState);
      }

      if (otherSearchState.closed) return;

      var tentativeDistance = cameFrom.g2 + distance(cameFrom.node, otherNode, link);

      if (tentativeDistance < otherSearchState.g2) {
        otherSearchState.g2 = tentativeDistance;
        otherSearchState.f2 = tentativeDistance + heuristic(from, otherSearchState.node);
        otherSearchState.p2 = cameFrom;
        if (otherSearchState.h2 < 0) {
          open2Set.push(otherSearchState);
        } else {
          open2Set.updateItem(otherSearchState.h2);
        }
      }
      var potentialMin = otherSearchState.g1 + otherSearchState.g2;
      if (potentialMin < lMin) {
        lMin = potentialMin;
        minNode = otherSearchState;
      }
    }

    function visitN2Oriented(otherNode, link) {
      // we are going backwards, graph needs to be reversed. 
      if (link.toId === cameFrom.node.id) return visitN2(otherNode, link);
    }
    function visitN1Oriented(otherNode, link) {
      // this is forward direction, so we should be coming FROM:
      if (link.fromId === cameFrom.node.id) return visitN1(otherNode, link);
    }
  }
}

function reconstructPath(searchState) {
  if (!searchState) return NO_PATH;

  var path = [searchState.node];
  var parent = searchState.p1;

  while (parent) {
    path.push(parent.node);
    parent = parent.p1;
  }

  var child = searchState.p2;

  while (child) {
    path.unshift(child.node);
    child = child.p2;
  }
  return path;
}

},{"../NodeHeap":56,"../defaultSettings.js":59,"../heuristics":60,"./makeNBASearchStatePool.js":63}],63:[function(require,module,exports){
module.exports = makeNBASearchStatePool;

/**
 * Creates new instance of NBASearchState. The instance stores information
 * about search state, and is used by NBA* algorithm.
 *
 * @param {Object} node - original graph node
 */
function NBASearchState(node) {
  /**
   * Original graph node.
   */
  this.node = node;

  /**
   * Parent of this node in forward search
   */
  this.p1 = null;

  /**
   * Parent of this node in reverse search
   */
  this.p2 = null;

  /**
   * If this is set to true, then the node was already processed
   * and we should not touch it anymore.
   */
  this.closed = false;

  /**
   * Actual distance from this node to its parent in forward search
   */
  this.g1 = Number.POSITIVE_INFINITY;

  /**
   * Actual distance from this node to its parent in reverse search
   */
  this.g2 = Number.POSITIVE_INFINITY;


  /**
   * Underestimated distance from this node to the path-finding source.
   */
  this.f1 = Number.POSITIVE_INFINITY;

  /**
   * Underestimated distance from this node to the path-finding target.
   */
  this.f2 = Number.POSITIVE_INFINITY;

  // used to reconstruct heap when fScore is updated. TODO: do I need them both?

  /**
   * Index of this node in the forward heap.
   */
  this.h1 = -1;

  /**
   * Index of this node in the reverse heap.
   */
  this.h2 = -1;
}

/**
 * As path-finding is memory-intensive process, we want to reduce pressure on
 * garbage collector. This class helps us to recycle path-finding nodes and significantly
 * reduces the search time (~20% faster than without it).
 */
function makeNBASearchStatePool() {
  var currentInCache = 0;
  var nodeCache = [];

  return {
    /**
     * Creates a new NBASearchState instance
     */
    createNewState: createNewState,

    /**
     * Marks all created instances available for recycling.
     */
    reset: reset
  };

  function reset() {
    currentInCache = 0;
  }

  function createNewState(node) {
    var cached = nodeCache[currentInCache];
    if (cached) {
      // TODO: This almost duplicates constructor code. Not sure if
      // it would impact performance if I move this code into a function
      cached.node = node;

      // How we came to this node?
      cached.p1 = null;
      cached.p2 = null;

      cached.closed = false;

      cached.g1 = Number.POSITIVE_INFINITY;
      cached.g2 = Number.POSITIVE_INFINITY;
      cached.f1 = Number.POSITIVE_INFINITY;
      cached.f2 = Number.POSITIVE_INFINITY;

      // used to reconstruct heap when fScore is updated.
      cached.h1 = -1;
      cached.h2 = -1;
    } else {
      cached = new NBASearchState(node);
      nodeCache[currentInCache] = cached;
    }
    currentInCache++;
    return cached;
  }
}

},{}],64:[function(require,module,exports){
module.exports = {
  aStar: require('./a-star/a-star.js'),
  aGreedy: require('./a-star/a-greedy-star'),
  nba: require('./a-star/nba/index.js'),
}

},{"./a-star/a-greedy-star":57,"./a-star/a-star.js":58,"./a-star/nba/index.js":62}],65:[function(require,module,exports){
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
	typeof define === 'function' && define.amd ? define(factory) :
	(global.quickselect = factory());
}(this, (function () { 'use strict';

function quickselect(arr, k, left, right, compare) {
    quickselectStep(arr, k, left || 0, right || (arr.length - 1), compare || defaultCompare);
}

function quickselectStep(arr, k, left, right, compare) {

    while (right > left) {
        if (right - left > 600) {
            var n = right - left + 1;
            var m = k - left + 1;
            var z = Math.log(n);
            var s = 0.5 * Math.exp(2 * z / 3);
            var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
            var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
            var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
            quickselectStep(arr, k, newLeft, newRight, compare);
        }

        var t = arr[k];
        var i = left;
        var j = right;

        swap(arr, left, k);
        if (compare(arr[right], t) > 0) swap(arr, left, right);

        while (i < j) {
            swap(arr, i, j);
            i++;
            j--;
            while (compare(arr[i], t) < 0) i++;
            while (compare(arr[j], t) > 0) j--;
        }

        if (compare(arr[left], t) === 0) swap(arr, left, j);
        else {
            j++;
            swap(arr, j, right);
        }

        if (j <= k) left = j + 1;
        if (k <= j) right = j - 1;
    }
}

function swap(arr, i, j) {
    var tmp = arr[i];
    arr[i] = arr[j];
    arr[j] = tmp;
}

function defaultCompare(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
}

return quickselect;

})));

},{}],66:[function(require,module,exports){
'use strict';

module.exports = rbush;
module.exports.default = rbush;

var quickselect = require('quickselect');

function rbush(maxEntries, format) {
    if (!(this instanceof rbush)) return new rbush(maxEntries, format);

    // max entries in a node is 9 by default; min node fill is 40% for best performance
    this._maxEntries = Math.max(4, maxEntries || 9);
    this._minEntries = Math.max(2, Math.ceil(this._maxEntries * 0.4));

    if (format) {
        this._initFormat(format);
    }

    this.clear();
}

rbush.prototype = {

    all: function () {
        return this._all(this.data, []);
    },

    search: function (bbox) {

        var node = this.data,
            result = [],
            toBBox = this.toBBox;

        if (!intersects(bbox, node)) return result;

        var nodesToSearch = [],
            i, len, child, childBBox;

        while (node) {
            for (i = 0, len = node.children.length; i < len; i++) {

                child = node.children[i];
                childBBox = node.leaf ? toBBox(child) : child;

                if (intersects(bbox, childBBox)) {
                    if (node.leaf) result.push(child);
                    else if (contains(bbox, childBBox)) this._all(child, result);
                    else nodesToSearch.push(child);
                }
            }
            node = nodesToSearch.pop();
        }

        return result;
    },

    collides: function (bbox) {

        var node = this.data,
            toBBox = this.toBBox;

        if (!intersects(bbox, node)) return false;

        var nodesToSearch = [],
            i, len, child, childBBox;

        while (node) {
            for (i = 0, len = node.children.length; i < len; i++) {

                child = node.children[i];
                childBBox = node.leaf ? toBBox(child) : child;

                if (intersects(bbox, childBBox)) {
                    if (node.leaf || contains(bbox, childBBox)) return true;
                    nodesToSearch.push(child);
                }
            }
            node = nodesToSearch.pop();
        }

        return false;
    },

    load: function (data) {
        if (!(data && data.length)) return this;

        if (data.length < this._minEntries) {
            for (var i = 0, len = data.length; i < len; i++) {
                this.insert(data[i]);
            }
            return this;
        }

        // recursively build the tree with the given data from scratch using OMT algorithm
        var node = this._build(data.slice(), 0, data.length - 1, 0);

        if (!this.data.children.length) {
            // save as is if tree is empty
            this.data = node;

        } else if (this.data.height === node.height) {
            // split root if trees have the same height
            this._splitRoot(this.data, node);

        } else {
            if (this.data.height < node.height) {
                // swap trees if inserted one is bigger
                var tmpNode = this.data;
                this.data = node;
                node = tmpNode;
            }

            // insert the small tree into the large tree at appropriate level
            this._insert(node, this.data.height - node.height - 1, true);
        }

        return this;
    },

    insert: function (item) {
        if (item) this._insert(item, this.data.height - 1);
        return this;
    },

    clear: function () {
        this.data = createNode([]);
        return this;
    },

    remove: function (item, equalsFn) {
        if (!item) return this;

        var node = this.data,
            bbox = this.toBBox(item),
            path = [],
            indexes = [],
            i, parent, index, goingUp;

        // depth-first iterative tree traversal
        while (node || path.length) {

            if (!node) { // go up
                node = path.pop();
                parent = path[path.length - 1];
                i = indexes.pop();
                goingUp = true;
            }

            if (node.leaf) { // check current node
                index = findItem(item, node.children, equalsFn);

                if (index !== -1) {
                    // item found, remove the item and condense tree upwards
                    node.children.splice(index, 1);
                    path.push(node);
                    this._condense(path);
                    return this;
                }
            }

            if (!goingUp && !node.leaf && contains(node, bbox)) { // go down
                path.push(node);
                indexes.push(i);
                i = 0;
                parent = node;
                node = node.children[0];

            } else if (parent) { // go right
                i++;
                node = parent.children[i];
                goingUp = false;

            } else node = null; // nothing found
        }

        return this;
    },

    toBBox: function (item) { return item; },

    compareMinX: compareNodeMinX,
    compareMinY: compareNodeMinY,

    toJSON: function () { return this.data; },

    fromJSON: function (data) {
        this.data = data;
        return this;
    },

    _all: function (node, result) {
        var nodesToSearch = [];
        while (node) {
            if (node.leaf) result.push.apply(result, node.children);
            else nodesToSearch.push.apply(nodesToSearch, node.children);

            node = nodesToSearch.pop();
        }
        return result;
    },

    _build: function (items, left, right, height) {

        var N = right - left + 1,
            M = this._maxEntries,
            node;

        if (N <= M) {
            // reached leaf level; return leaf
            node = createNode(items.slice(left, right + 1));
            calcBBox(node, this.toBBox);
            return node;
        }

        if (!height) {
            // target height of the bulk-loaded tree
            height = Math.ceil(Math.log(N) / Math.log(M));

            // target number of root entries to maximize storage utilization
            M = Math.ceil(N / Math.pow(M, height - 1));
        }

        node = createNode([]);
        node.leaf = false;
        node.height = height;

        // split the items into M mostly square tiles

        var N2 = Math.ceil(N / M),
            N1 = N2 * Math.ceil(Math.sqrt(M)),
            i, j, right2, right3;

        multiSelect(items, left, right, N1, this.compareMinX);

        for (i = left; i <= right; i += N1) {

            right2 = Math.min(i + N1 - 1, right);

            multiSelect(items, i, right2, N2, this.compareMinY);

            for (j = i; j <= right2; j += N2) {

                right3 = Math.min(j + N2 - 1, right2);

                // pack each entry recursively
                node.children.push(this._build(items, j, right3, height - 1));
            }
        }

        calcBBox(node, this.toBBox);

        return node;
    },

    _chooseSubtree: function (bbox, node, level, path) {

        var i, len, child, targetNode, area, enlargement, minArea, minEnlargement;

        while (true) {
            path.push(node);

            if (node.leaf || path.length - 1 === level) break;

            minArea = minEnlargement = Infinity;

            for (i = 0, len = node.children.length; i < len; i++) {
                child = node.children[i];
                area = bboxArea(child);
                enlargement = enlargedArea(bbox, child) - area;

                // choose entry with the least area enlargement
                if (enlargement < minEnlargement) {
                    minEnlargement = enlargement;
                    minArea = area < minArea ? area : minArea;
                    targetNode = child;

                } else if (enlargement === minEnlargement) {
                    // otherwise choose one with the smallest area
                    if (area < minArea) {
                        minArea = area;
                        targetNode = child;
                    }
                }
            }

            node = targetNode || node.children[0];
        }

        return node;
    },

    _insert: function (item, level, isNode) {

        var toBBox = this.toBBox,
            bbox = isNode ? item : toBBox(item),
            insertPath = [];

        // find the best node for accommodating the item, saving all nodes along the path too
        var node = this._chooseSubtree(bbox, this.data, level, insertPath);

        // put the item into the node
        node.children.push(item);
        extend(node, bbox);

        // split on node overflow; propagate upwards if necessary
        while (level >= 0) {
            if (insertPath[level].children.length > this._maxEntries) {
                this._split(insertPath, level);
                level--;
            } else break;
        }

        // adjust bboxes along the insertion path
        this._adjustParentBBoxes(bbox, insertPath, level);
    },

    // split overflowed node into two
    _split: function (insertPath, level) {

        var node = insertPath[level],
            M = node.children.length,
            m = this._minEntries;

        this._chooseSplitAxis(node, m, M);

        var splitIndex = this._chooseSplitIndex(node, m, M);

        var newNode = createNode(node.children.splice(splitIndex, node.children.length - splitIndex));
        newNode.height = node.height;
        newNode.leaf = node.leaf;

        calcBBox(node, this.toBBox);
        calcBBox(newNode, this.toBBox);

        if (level) insertPath[level - 1].children.push(newNode);
        else this._splitRoot(node, newNode);
    },

    _splitRoot: function (node, newNode) {
        // split root node
        this.data = createNode([node, newNode]);
        this.data.height = node.height + 1;
        this.data.leaf = false;
        calcBBox(this.data, this.toBBox);
    },

    _chooseSplitIndex: function (node, m, M) {

        var i, bbox1, bbox2, overlap, area, minOverlap, minArea, index;

        minOverlap = minArea = Infinity;

        for (i = m; i <= M - m; i++) {
            bbox1 = distBBox(node, 0, i, this.toBBox);
            bbox2 = distBBox(node, i, M, this.toBBox);

            overlap = intersectionArea(bbox1, bbox2);
            area = bboxArea(bbox1) + bboxArea(bbox2);

            // choose distribution with minimum overlap
            if (overlap < minOverlap) {
                minOverlap = overlap;
                index = i;

                minArea = area < minArea ? area : minArea;

            } else if (overlap === minOverlap) {
                // otherwise choose distribution with minimum area
                if (area < minArea) {
                    minArea = area;
                    index = i;
                }
            }
        }

        return index;
    },

    // sorts node children by the best axis for split
    _chooseSplitAxis: function (node, m, M) {

        var compareMinX = node.leaf ? this.compareMinX : compareNodeMinX,
            compareMinY = node.leaf ? this.compareMinY : compareNodeMinY,
            xMargin = this._allDistMargin(node, m, M, compareMinX),
            yMargin = this._allDistMargin(node, m, M, compareMinY);

        // if total distributions margin value is minimal for x, sort by minX,
        // otherwise it's already sorted by minY
        if (xMargin < yMargin) node.children.sort(compareMinX);
    },

    // total margin of all possible split distributions where each node is at least m full
    _allDistMargin: function (node, m, M, compare) {

        node.children.sort(compare);

        var toBBox = this.toBBox,
            leftBBox = distBBox(node, 0, m, toBBox),
            rightBBox = distBBox(node, M - m, M, toBBox),
            margin = bboxMargin(leftBBox) + bboxMargin(rightBBox),
            i, child;

        for (i = m; i < M - m; i++) {
            child = node.children[i];
            extend(leftBBox, node.leaf ? toBBox(child) : child);
            margin += bboxMargin(leftBBox);
        }

        for (i = M - m - 1; i >= m; i--) {
            child = node.children[i];
            extend(rightBBox, node.leaf ? toBBox(child) : child);
            margin += bboxMargin(rightBBox);
        }

        return margin;
    },

    _adjustParentBBoxes: function (bbox, path, level) {
        // adjust bboxes along the given tree path
        for (var i = level; i >= 0; i--) {
            extend(path[i], bbox);
        }
    },

    _condense: function (path) {
        // go through the path, removing empty nodes and updating bboxes
        for (var i = path.length - 1, siblings; i >= 0; i--) {
            if (path[i].children.length === 0) {
                if (i > 0) {
                    siblings = path[i - 1].children;
                    siblings.splice(siblings.indexOf(path[i]), 1);

                } else this.clear();

            } else calcBBox(path[i], this.toBBox);
        }
    },

    _initFormat: function (format) {
        // data format (minX, minY, maxX, maxY accessors)

        // uses eval-type function compilation instead of just accepting a toBBox function
        // because the algorithms are very sensitive to sorting functions performance,
        // so they should be dead simple and without inner calls

        var compareArr = ['return a', ' - b', ';'];

        this.compareMinX = new Function('a', 'b', compareArr.join(format[0]));
        this.compareMinY = new Function('a', 'b', compareArr.join(format[1]));

        this.toBBox = new Function('a',
            'return {minX: a' + format[0] +
            ', minY: a' + format[1] +
            ', maxX: a' + format[2] +
            ', maxY: a' + format[3] + '};');
    }
};

function findItem(item, items, equalsFn) {
    if (!equalsFn) return items.indexOf(item);

    for (var i = 0; i < items.length; i++) {
        if (equalsFn(item, items[i])) return i;
    }
    return -1;
}

// calculate node's bbox from bboxes of its children
function calcBBox(node, toBBox) {
    distBBox(node, 0, node.children.length, toBBox, node);
}

// min bounding rectangle of node children from k to p-1
function distBBox(node, k, p, toBBox, destNode) {
    if (!destNode) destNode = createNode(null);
    destNode.minX = Infinity;
    destNode.minY = Infinity;
    destNode.maxX = -Infinity;
    destNode.maxY = -Infinity;

    for (var i = k, child; i < p; i++) {
        child = node.children[i];
        extend(destNode, node.leaf ? toBBox(child) : child);
    }

    return destNode;
}

function extend(a, b) {
    a.minX = Math.min(a.minX, b.minX);
    a.minY = Math.min(a.minY, b.minY);
    a.maxX = Math.max(a.maxX, b.maxX);
    a.maxY = Math.max(a.maxY, b.maxY);
    return a;
}

function compareNodeMinX(a, b) { return a.minX - b.minX; }
function compareNodeMinY(a, b) { return a.minY - b.minY; }

function bboxArea(a)   { return (a.maxX - a.minX) * (a.maxY - a.minY); }
function bboxMargin(a) { return (a.maxX - a.minX) + (a.maxY - a.minY); }

function enlargedArea(a, b) {
    return (Math.max(b.maxX, a.maxX) - Math.min(b.minX, a.minX)) *
           (Math.max(b.maxY, a.maxY) - Math.min(b.minY, a.minY));
}

function intersectionArea(a, b) {
    var minX = Math.max(a.minX, b.minX),
        minY = Math.max(a.minY, b.minY),
        maxX = Math.min(a.maxX, b.maxX),
        maxY = Math.min(a.maxY, b.maxY);

    return Math.max(0, maxX - minX) *
           Math.max(0, maxY - minY);
}

function contains(a, b) {
    return a.minX <= b.minX &&
           a.minY <= b.minY &&
           b.maxX <= a.maxX &&
           b.maxY <= a.maxY;
}

function intersects(a, b) {
    return b.minX <= a.maxX &&
           b.minY <= a.maxY &&
           b.maxX >= a.minX &&
           b.maxY >= a.minY;
}

function createNode(children) {
    return {
        children: children,
        height: 1,
        leaf: true,
        minX: Infinity,
        minY: Infinity,
        maxX: -Infinity,
        maxY: -Infinity
    };
}

// sort an array so that items come in groups of n unsorted items, with groups sorted between each other;
// combines selection algorithm with binary divide & conquer approach

function multiSelect(arr, left, right, n, compare) {
    var stack = [left, right],
        mid;

    while (stack.length) {
        right = stack.pop();
        left = stack.pop();

        if (right - left <= n) continue;

        mid = left + Math.ceil((right - left) / n / 2) * n;
        quickselect(arr, mid, left, right, compare);

        stack.push(left, mid, mid, right);
    }
}

},{"quickselect":65}],67:[function(require,module,exports){
!function(t,e){"object"==typeof exports&&"undefined"!=typeof module?e(exports):"function"==typeof define&&define.amd?define(["exports"],e):e(t.jsts={})}(this,function(t){"use strict";function e(){}function n(t){this.message=t||""}function i(t){this.message=t||""}function r(t){this.message=t||""}function o(){}function s(t){return null===t?Mt:t.color}function a(t){return null===t?null:t.parent}function u(t,e){null!==t&&(t.color=e)}function l(t){return null===t?null:t.left}function c(t){return null===t?null:t.right}function p(){this.root_=null,this.size_=0}function h(){}function f(){this.array_=[],arguments[0]instanceof It&&this.addAll(arguments[0])}function g(){}function d(t){this.message=t||""}function y(){this.array_=[]}"fill"in Array.prototype||Object.defineProperty(Array.prototype,"fill",{configurable:!0,value:function(t){if(void 0===this||null===this)throw new TypeError(this+" is not an object");var e=Object(this),n=Math.max(Math.min(e.length,9007199254740991),0)||0,i=1 in arguments?parseInt(Number(arguments[1]),10)||0:0;i=i<0?Math.max(n+i,0):Math.min(i,n);var r=2 in arguments&&void 0!==arguments[2]?parseInt(Number(arguments[2]),10)||0:n;for(r=r<0?Math.max(n+arguments[2],0):Math.min(r,n);i<r;)e[i]=t,++i;return e},writable:!0}),Number.isFinite=Number.isFinite||function(t){return"number"==typeof t&&isFinite(t)},Number.isInteger=Number.isInteger||function(t){return"number"==typeof t&&isFinite(t)&&Math.floor(t)===t},Number.parseFloat=Number.parseFloat||parseFloat,Number.isNaN=Number.isNaN||function(t){return t!=t},Math.trunc=Math.trunc||function(t){return t<0?Math.ceil(t):Math.floor(t)};var _=function(){};_.prototype.interfaces_=function(){return[]},_.prototype.getClass=function(){return _},_.prototype.equalsWithTolerance=function(t,e,n){return Math.abs(t-e)<=n};var m=function(t){function e(e){t.call(this,e),this.name="IllegalArgumentException",this.message=e,this.stack=(new t).stack}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e}(Error),v=function(){},I={MAX_VALUE:{configurable:!0}};v.isNaN=function(t){return Number.isNaN(t)},v.doubleToLongBits=function(t){return t},v.longBitsToDouble=function(t){return t},v.isInfinite=function(t){return!Number.isFinite(t)},I.MAX_VALUE.get=function(){return Number.MAX_VALUE},Object.defineProperties(v,I);var E=function(){},x=function(){},N=function(){},C=function t(){if(this.x=null,this.y=null,this.z=null,0===arguments.length)this.x=0,this.y=0,this.z=t.NULL_ORDINATE;else if(1===arguments.length){var e=arguments[0];this.x=e.x,this.y=e.y,this.z=e.z}else 2===arguments.length?(this.x=arguments[0],this.y=arguments[1],this.z=t.NULL_ORDINATE):3===arguments.length&&(this.x=arguments[0],this.y=arguments[1],this.z=arguments[2])},S={DimensionalComparator:{configurable:!0},serialVersionUID:{configurable:!0},NULL_ORDINATE:{configurable:!0},X:{configurable:!0},Y:{configurable:!0},Z:{configurable:!0}};C.prototype.setOrdinate=function(t,e){switch(t){case C.X:this.x=e;break;case C.Y:this.y=e;break;case C.Z:this.z=e;break;default:throw new m("Invalid ordinate index: "+t)}},C.prototype.equals2D=function(){if(1===arguments.length){var t=arguments[0];return this.x===t.x&&this.y===t.y}if(2===arguments.length){var e=arguments[0],n=arguments[1];return!!_.equalsWithTolerance(this.x,e.x,n)&&!!_.equalsWithTolerance(this.y,e.y,n)}},C.prototype.getOrdinate=function(t){switch(t){case C.X:return this.x;case C.Y:return this.y;case C.Z:return this.z}throw new m("Invalid ordinate index: "+t)},C.prototype.equals3D=function(t){return this.x===t.x&&this.y===t.y&&(this.z===t.z||v.isNaN(this.z))&&v.isNaN(t.z)},C.prototype.equals=function(t){return t instanceof C&&this.equals2D(t)},C.prototype.equalInZ=function(t,e){return _.equalsWithTolerance(this.z,t.z,e)},C.prototype.compareTo=function(t){var e=t;return this.x<e.x?-1:this.x>e.x?1:this.y<e.y?-1:this.y>e.y?1:0},C.prototype.clone=function(){},C.prototype.copy=function(){return new C(this)},C.prototype.toString=function(){return"("+this.x+", "+this.y+", "+this.z+")"},C.prototype.distance3D=function(t){var e=this.x-t.x,n=this.y-t.y,i=this.z-t.z;return Math.sqrt(e*e+n*n+i*i)},C.prototype.distance=function(t){var e=this.x-t.x,n=this.y-t.y;return Math.sqrt(e*e+n*n)},C.prototype.hashCode=function(){var t=17;return t=37*t+C.hashCode(this.x),t=37*t+C.hashCode(this.y)},C.prototype.setCoordinate=function(t){this.x=t.x,this.y=t.y,this.z=t.z},C.prototype.interfaces_=function(){return[E,x,e]},C.prototype.getClass=function(){return C},C.hashCode=function(){if(1===arguments.length){var t=arguments[0],e=v.doubleToLongBits(t);return Math.trunc((e^e)>>>32)}},S.DimensionalComparator.get=function(){return L},S.serialVersionUID.get=function(){return 0x5cbf2c235c7e5800},S.NULL_ORDINATE.get=function(){return v.NaN},S.X.get=function(){return 0},S.Y.get=function(){return 1},S.Z.get=function(){return 2},Object.defineProperties(C,S);var L=function(t){if(this._dimensionsToTest=2,0===arguments.length);else if(1===arguments.length){var e=arguments[0];if(2!==e&&3!==e)throw new m("only 2 or 3 dimensions may be specified");this._dimensionsToTest=e}};L.prototype.compare=function(t,e){var n=t,i=e,r=L.compare(n.x,i.x);if(0!==r)return r;var o=L.compare(n.y,i.y);if(0!==o)return o;if(this._dimensionsToTest<=2)return 0;return L.compare(n.z,i.z)},L.prototype.interfaces_=function(){return[N]},L.prototype.getClass=function(){return L},L.compare=function(t,e){return t<e?-1:t>e?1:v.isNaN(t)?v.isNaN(e)?0:-1:v.isNaN(e)?1:0};var b=function(){};b.prototype.create=function(){},b.prototype.interfaces_=function(){return[]},b.prototype.getClass=function(){return b};var w=function(){},O={INTERIOR:{configurable:!0},BOUNDARY:{configurable:!0},EXTERIOR:{configurable:!0},NONE:{configurable:!0}};w.prototype.interfaces_=function(){return[]},w.prototype.getClass=function(){return w},w.toLocationSymbol=function(t){switch(t){case w.EXTERIOR:return"e";case w.BOUNDARY:return"b";case w.INTERIOR:return"i";case w.NONE:return"-"}throw new m("Unknown location value: "+t)},O.INTERIOR.get=function(){return 0},O.BOUNDARY.get=function(){return 1},O.EXTERIOR.get=function(){return 2},O.NONE.get=function(){return-1},Object.defineProperties(w,O);var T=function(t,e){return t.interfaces_&&t.interfaces_().indexOf(e)>-1},R=function(){},P={LOG_10:{configurable:!0}};R.prototype.interfaces_=function(){return[]},R.prototype.getClass=function(){return R},R.log10=function(t){var e=Math.log(t);return v.isInfinite(e)?e:v.isNaN(e)?e:e/R.LOG_10},R.min=function(t,e,n,i){var r=t;return e<r&&(r=e),n<r&&(r=n),i<r&&(r=i),r},R.clamp=function(){if("number"==typeof arguments[2]&&"number"==typeof arguments[0]&&"number"==typeof arguments[1]){var t=arguments[0],e=arguments[1],n=arguments[2];return t<e?e:t>n?n:t}if(Number.isInteger(arguments[2])&&Number.isInteger(arguments[0])&&Number.isInteger(arguments[1])){var i=arguments[0],r=arguments[1],o=arguments[2];return i<r?r:i>o?o:i}},R.wrap=function(t,e){return t<0?e- -t%e:t%e},R.max=function(){if(3===arguments.length){var t=arguments[0],e=arguments[1],n=arguments[2],i=t;return e>i&&(i=e),n>i&&(i=n),i}if(4===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2],a=arguments[3],u=r;return o>u&&(u=o),s>u&&(u=s),a>u&&(u=a),u}},R.average=function(t,e){return(t+e)/2},P.LOG_10.get=function(){return Math.log(10)},Object.defineProperties(R,P);var D=function(t){this.str=t};D.prototype.append=function(t){this.str+=t},D.prototype.setCharAt=function(t,e){this.str=this.str.substr(0,t)+e+this.str.substr(t+1)},D.prototype.toString=function(t){return this.str};var M=function(t){this.value=t};M.prototype.intValue=function(){return this.value},M.prototype.compareTo=function(t){return this.value<t?-1:this.value>t?1:0},M.isNaN=function(t){return Number.isNaN(t)};var A=function(){};A.isWhitespace=function(t){return t<=32&&t>=0||127===t},A.toUpperCase=function(t){return t.toUpperCase()};var F=function t(){if(this._hi=0,this._lo=0,0===arguments.length)this.init(0);else if(1===arguments.length){if("number"==typeof arguments[0]){var e=arguments[0];this.init(e)}else if(arguments[0]instanceof t){var n=arguments[0];this.init(n)}else if("string"==typeof arguments[0]){var i=arguments[0];t.call(this,t.parse(i))}}else if(2===arguments.length){var r=arguments[0],o=arguments[1];this.init(r,o)}},G={PI:{configurable:!0},TWO_PI:{configurable:!0},PI_2:{configurable:!0},E:{configurable:!0},NaN:{configurable:!0},EPS:{configurable:!0},SPLIT:{configurable:!0},MAX_PRINT_DIGITS:{configurable:!0},TEN:{configurable:!0},ONE:{configurable:!0},SCI_NOT_EXPONENT_CHAR:{configurable:!0},SCI_NOT_ZERO:{configurable:!0}};F.prototype.le=function(t){return(this._hi<t._hi||this._hi===t._hi)&&this._lo<=t._lo},F.prototype.extractSignificantDigits=function(t,e){var n=this.abs(),i=F.magnitude(n._hi),r=F.TEN.pow(i);(n=n.divide(r)).gt(F.TEN)?(n=n.divide(F.TEN),i+=1):n.lt(F.ONE)&&(n=n.multiply(F.TEN),i-=1);for(var o=i+1,s=new D,a=F.MAX_PRINT_DIGITS-1,u=0;u<=a;u++){t&&u===o&&s.append(".");var l=Math.trunc(n._hi);if(l<0)break;var c=!1,p=0;l>9?(c=!0,p="9"):p="0"+l,s.append(p),n=n.subtract(F.valueOf(l)).multiply(F.TEN),c&&n.selfAdd(F.TEN);var h=!0,f=F.magnitude(n._hi);if(f<0&&Math.abs(f)>=a-u&&(h=!1),!h)break}return e[0]=i,s.toString()},F.prototype.sqr=function(){return this.multiply(this)},F.prototype.doubleValue=function(){return this._hi+this._lo},F.prototype.subtract=function(){if(arguments[0]instanceof F){var t=arguments[0];return this.add(t.negate())}if("number"==typeof arguments[0]){var e=arguments[0];return this.add(-e)}},F.prototype.equals=function(){if(1===arguments.length){var t=arguments[0];return this._hi===t._hi&&this._lo===t._lo}},F.prototype.isZero=function(){return 0===this._hi&&0===this._lo},F.prototype.selfSubtract=function(){if(arguments[0]instanceof F){var t=arguments[0];return this.isNaN()?this:this.selfAdd(-t._hi,-t._lo)}if("number"==typeof arguments[0]){var e=arguments[0];return this.isNaN()?this:this.selfAdd(-e,0)}},F.prototype.getSpecialNumberString=function(){return this.isZero()?"0.0":this.isNaN()?"NaN ":null},F.prototype.min=function(t){return this.le(t)?this:t},F.prototype.selfDivide=function(){if(1===arguments.length){if(arguments[0]instanceof F){var t=arguments[0];return this.selfDivide(t._hi,t._lo)}if("number"==typeof arguments[0]){var e=arguments[0];return this.selfDivide(e,0)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1],r=null,o=null,s=null,a=null,u=null,l=null,c=null,p=null;return u=this._hi/n,l=F.SPLIT*u,r=l-u,p=F.SPLIT*n,r=l-r,o=u-r,s=p-n,c=u*n,s=p-s,a=n-s,p=r*s-c+r*a+o*s+o*a,l=(this._hi-c-p+this._lo-u*i)/n,p=u+l,this._hi=p,this._lo=u-p+l,this}},F.prototype.dump=function(){return"DD<"+this._hi+", "+this._lo+">"},F.prototype.divide=function(){if(arguments[0]instanceof F){var t=arguments[0],e=null,n=null,i=null,r=null,o=null,s=null,a=null,u=null;n=(o=this._hi/t._hi)-(e=(s=F.SPLIT*o)-(e=s-o)),u=e*(i=(u=F.SPLIT*t._hi)-(i=u-t._hi))-(a=o*t._hi)+e*(r=t._hi-i)+n*i+n*r,s=(this._hi-a-u+this._lo-o*t._lo)/t._hi;return new F(u=o+s,o-u+s)}if("number"==typeof arguments[0]){var l=arguments[0];return v.isNaN(l)?F.createNaN():F.copy(this).selfDivide(l,0)}},F.prototype.ge=function(t){return(this._hi>t._hi||this._hi===t._hi)&&this._lo>=t._lo},F.prototype.pow=function(t){if(0===t)return F.valueOf(1);var e=new F(this),n=F.valueOf(1),i=Math.abs(t);if(i>1)for(;i>0;)i%2==1&&n.selfMultiply(e),(i/=2)>0&&(e=e.sqr());else n=e;return t<0?n.reciprocal():n},F.prototype.ceil=function(){if(this.isNaN())return F.NaN;var t=Math.ceil(this._hi),e=0;return t===this._hi&&(e=Math.ceil(this._lo)),new F(t,e)},F.prototype.compareTo=function(t){var e=t;return this._hi<e._hi?-1:this._hi>e._hi?1:this._lo<e._lo?-1:this._lo>e._lo?1:0},F.prototype.rint=function(){if(this.isNaN())return this;return this.add(.5).floor()},F.prototype.setValue=function(){if(arguments[0]instanceof F){var t=arguments[0];return this.init(t),this}if("number"==typeof arguments[0]){var e=arguments[0];return this.init(e),this}},F.prototype.max=function(t){return this.ge(t)?this:t},F.prototype.sqrt=function(){if(this.isZero())return F.valueOf(0);if(this.isNegative())return F.NaN;var t=1/Math.sqrt(this._hi),e=this._hi*t,n=F.valueOf(e),i=this.subtract(n.sqr())._hi*(.5*t);return n.add(i)},F.prototype.selfAdd=function(){if(1===arguments.length){if(arguments[0]instanceof F){var t=arguments[0];return this.selfAdd(t._hi,t._lo)}if("number"==typeof arguments[0]){var e=arguments[0],n=null,i=null,r=null,o=null,s=null,a=null;return r=this._hi+e,s=r-this._hi,o=r-s,o=e-s+(this._hi-o),a=o+this._lo,n=r+a,i=a+(r-n),this._hi=n+i,this._lo=i+(n-this._hi),this}}else if(2===arguments.length){var u=arguments[0],l=arguments[1],c=null,p=null,h=null,f=null,g=null,d=null,y=null;f=this._hi+u,p=this._lo+l,g=f-(d=f-this._hi),h=p-(y=p-this._lo);var _=(c=f+(d=(g=u-d+(this._hi-g))+p))+(d=(h=l-y+(this._lo-h))+(d+(f-c))),m=d+(c-_);return this._hi=_,this._lo=m,this}},F.prototype.selfMultiply=function(){if(1===arguments.length){if(arguments[0]instanceof F){var t=arguments[0];return this.selfMultiply(t._hi,t._lo)}if("number"==typeof arguments[0]){var e=arguments[0];return this.selfMultiply(e,0)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1],r=null,o=null,s=null,a=null,u=null,l=null;r=(u=F.SPLIT*this._hi)-this._hi,l=F.SPLIT*n,r=u-r,o=this._hi-r,s=l-n;var c=(u=this._hi*n)+(l=r*(s=l-s)-u+r*(a=n-s)+o*s+o*a+(this._hi*i+this._lo*n)),p=l+(r=u-c);return this._hi=c,this._lo=p,this}},F.prototype.selfSqr=function(){return this.selfMultiply(this)},F.prototype.floor=function(){if(this.isNaN())return F.NaN;var t=Math.floor(this._hi),e=0;return t===this._hi&&(e=Math.floor(this._lo)),new F(t,e)},F.prototype.negate=function(){return this.isNaN()?this:new F(-this._hi,-this._lo)},F.prototype.clone=function(){},F.prototype.multiply=function(){if(arguments[0]instanceof F){var t=arguments[0];return t.isNaN()?F.createNaN():F.copy(this).selfMultiply(t)}if("number"==typeof arguments[0]){var e=arguments[0];return v.isNaN(e)?F.createNaN():F.copy(this).selfMultiply(e,0)}},F.prototype.isNaN=function(){return v.isNaN(this._hi)},F.prototype.intValue=function(){return Math.trunc(this._hi)},F.prototype.toString=function(){var t=F.magnitude(this._hi);return t>=-3&&t<=20?this.toStandardNotation():this.toSciNotation()},F.prototype.toStandardNotation=function(){var t=this.getSpecialNumberString();if(null!==t)return t;var e=new Array(1).fill(null),n=this.extractSignificantDigits(!0,e),i=e[0]+1,r=n;if("."===n.charAt(0))r="0"+n;else if(i<0)r="0."+F.stringOfChar("0",-i)+n;else if(-1===n.indexOf(".")){var o=i-n.length;r=n+F.stringOfChar("0",o)+".0"}return this.isNegative()?"-"+r:r},F.prototype.reciprocal=function(){var t=null,e=null,n=null,i=null,r=null,o=null,s=null,a=null;e=(r=1/this._hi)-(t=(o=F.SPLIT*r)-(t=o-r)),n=(a=F.SPLIT*this._hi)-this._hi;var u=r+(o=(1-(s=r*this._hi)-(a=t*(n=a-n)-s+t*(i=this._hi-n)+e*n+e*i)-r*this._lo)/this._hi);return new F(u,r-u+o)},F.prototype.toSciNotation=function(){if(this.isZero())return F.SCI_NOT_ZERO;var t=this.getSpecialNumberString();if(null!==t)return t;var e=new Array(1).fill(null),n=this.extractSignificantDigits(!1,e),i=F.SCI_NOT_EXPONENT_CHAR+e[0];if("0"===n.charAt(0))throw new Error("Found leading zero: "+n);var r="";n.length>1&&(r=n.substring(1));var o=n.charAt(0)+"."+r;return this.isNegative()?"-"+o+i:o+i},F.prototype.abs=function(){return this.isNaN()?F.NaN:this.isNegative()?this.negate():new F(this)},F.prototype.isPositive=function(){return(this._hi>0||0===this._hi)&&this._lo>0},F.prototype.lt=function(t){return(this._hi<t._hi||this._hi===t._hi)&&this._lo<t._lo},F.prototype.add=function(){if(arguments[0]instanceof F){var t=arguments[0];return F.copy(this).selfAdd(t)}if("number"==typeof arguments[0]){var e=arguments[0];return F.copy(this).selfAdd(e)}},F.prototype.init=function(){if(1===arguments.length){if("number"==typeof arguments[0]){var t=arguments[0];this._hi=t,this._lo=0}else if(arguments[0]instanceof F){var e=arguments[0];this._hi=e._hi,this._lo=e._lo}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this._hi=n,this._lo=i}},F.prototype.gt=function(t){return(this._hi>t._hi||this._hi===t._hi)&&this._lo>t._lo},F.prototype.isNegative=function(){return(this._hi<0||0===this._hi)&&this._lo<0},F.prototype.trunc=function(){return this.isNaN()?F.NaN:this.isPositive()?this.floor():this.ceil()},F.prototype.signum=function(){return this._hi>0?1:this._hi<0?-1:this._lo>0?1:this._lo<0?-1:0},F.prototype.interfaces_=function(){return[e,E,x]},F.prototype.getClass=function(){return F},F.sqr=function(t){return F.valueOf(t).selfMultiply(t)},F.valueOf=function(){if("string"==typeof arguments[0]){var t=arguments[0];return F.parse(t)}if("number"==typeof arguments[0]){var e=arguments[0];return new F(e)}},F.sqrt=function(t){return F.valueOf(t).sqrt()},F.parse=function(t){for(var e=0,n=t.length;A.isWhitespace(t.charAt(e));)e++;var i=!1;if(e<n){var r=t.charAt(e);"-"!==r&&"+"!==r||(e++,"-"===r&&(i=!0))}for(var o=new F,s=0,a=0,u=0;!(e>=n);){var l=t.charAt(e);if(e++,A.isDigit(l)){var c=l-"0";o.selfMultiply(F.TEN),o.selfAdd(c),s++}else{if("."!==l){if("e"===l||"E"===l){var p=t.substring(e);try{u=M.parseInt(p)}catch(e){throw e instanceof Error?new Error("Invalid exponent "+p+" in string "+t):e}break}throw new Error("Unexpected character '"+l+"' at position "+e+" in string "+t)}a=s}}var h=o,f=s-a-u;if(0===f)h=o;else if(f>0){var g=F.TEN.pow(f);h=o.divide(g)}else if(f<0){var d=F.TEN.pow(-f);h=o.multiply(d)}return i?h.negate():h},F.createNaN=function(){return new F(v.NaN,v.NaN)},F.copy=function(t){return new F(t)},F.magnitude=function(t){var e=Math.abs(t),n=Math.log(e)/Math.log(10),i=Math.trunc(Math.floor(n));return 10*Math.pow(10,i)<=e&&(i+=1),i},F.stringOfChar=function(t,e){for(var n=new D,i=0;i<e;i++)n.append(t);return n.toString()},G.PI.get=function(){return new F(3.141592653589793,1.2246467991473532e-16)},G.TWO_PI.get=function(){return new F(6.283185307179586,2.4492935982947064e-16)},G.PI_2.get=function(){return new F(1.5707963267948966,6.123233995736766e-17)},G.E.get=function(){return new F(2.718281828459045,1.4456468917292502e-16)},G.NaN.get=function(){return new F(v.NaN,v.NaN)},G.EPS.get=function(){return 1.23259516440783e-32},G.SPLIT.get=function(){return 134217729},G.MAX_PRINT_DIGITS.get=function(){return 32},G.TEN.get=function(){return F.valueOf(10)},G.ONE.get=function(){return F.valueOf(1)},G.SCI_NOT_EXPONENT_CHAR.get=function(){return"E"},G.SCI_NOT_ZERO.get=function(){return"0.0E0"},Object.defineProperties(F,G);var q=function(){},B={DP_SAFE_EPSILON:{configurable:!0}};q.prototype.interfaces_=function(){return[]},q.prototype.getClass=function(){return q},q.orientationIndex=function(t,e,n){var i=q.orientationIndexFilter(t,e,n);if(i<=1)return i;var r=F.valueOf(e.x).selfAdd(-t.x),o=F.valueOf(e.y).selfAdd(-t.y),s=F.valueOf(n.x).selfAdd(-e.x),a=F.valueOf(n.y).selfAdd(-e.y);return r.selfMultiply(a).selfSubtract(o.selfMultiply(s)).signum()},q.signOfDet2x2=function(t,e,n,i){return t.multiply(i).selfSubtract(e.multiply(n)).signum()},q.intersection=function(t,e,n,i){var r=F.valueOf(i.y).selfSubtract(n.y).selfMultiply(F.valueOf(e.x).selfSubtract(t.x)),o=F.valueOf(i.x).selfSubtract(n.x).selfMultiply(F.valueOf(e.y).selfSubtract(t.y)),s=r.subtract(o),a=F.valueOf(i.x).selfSubtract(n.x).selfMultiply(F.valueOf(t.y).selfSubtract(n.y)),u=F.valueOf(i.y).selfSubtract(n.y).selfMultiply(F.valueOf(t.x).selfSubtract(n.x)),l=a.subtract(u).selfDivide(s).doubleValue(),c=F.valueOf(t.x).selfAdd(F.valueOf(e.x).selfSubtract(t.x).selfMultiply(l)).doubleValue(),p=F.valueOf(e.x).selfSubtract(t.x).selfMultiply(F.valueOf(t.y).selfSubtract(n.y)),h=F.valueOf(e.y).selfSubtract(t.y).selfMultiply(F.valueOf(t.x).selfSubtract(n.x)),f=p.subtract(h).selfDivide(s).doubleValue(),g=F.valueOf(n.y).selfAdd(F.valueOf(i.y).selfSubtract(n.y).selfMultiply(f)).doubleValue();return new C(c,g)},q.orientationIndexFilter=function(t,e,n){var i=null,r=(t.x-n.x)*(e.y-n.y),o=(t.y-n.y)*(e.x-n.x),s=r-o;if(r>0){if(o<=0)return q.signum(s);i=r+o}else{if(!(r<0))return q.signum(s);if(o>=0)return q.signum(s);i=-r-o}var a=q.DP_SAFE_EPSILON*i;return s>=a||-s>=a?q.signum(s):2},q.signum=function(t){return t>0?1:t<0?-1:0},B.DP_SAFE_EPSILON.get=function(){return 1e-15},Object.defineProperties(q,B);var V=function(){},U={X:{configurable:!0},Y:{configurable:!0},Z:{configurable:!0},M:{configurable:!0}};U.X.get=function(){return 0},U.Y.get=function(){return 1},U.Z.get=function(){return 2},U.M.get=function(){return 3},V.prototype.setOrdinate=function(t,e,n){},V.prototype.size=function(){},V.prototype.getOrdinate=function(t,e){},V.prototype.getCoordinate=function(){},V.prototype.getCoordinateCopy=function(t){},V.prototype.getDimension=function(){},V.prototype.getX=function(t){},V.prototype.clone=function(){},V.prototype.expandEnvelope=function(t){},V.prototype.copy=function(){},V.prototype.getY=function(t){},V.prototype.toCoordinateArray=function(){},V.prototype.interfaces_=function(){return[x]},V.prototype.getClass=function(){return V},Object.defineProperties(V,U);var z=function(){},X=function(t){function e(){t.call(this,"Projective point not representable on the Cartesian plane.")}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(z),Y=function(){};Y.arraycopy=function(t,e,n,i,r){for(var o=0,s=e;s<e+r;s++)n[i+o]=t[s],o++},Y.getProperty=function(t){return{"line.separator":"\n"}[t]};var k=function t(){if(this.x=null,this.y=null,this.w=null,0===arguments.length)this.x=0,this.y=0,this.w=1;else if(1===arguments.length){var e=arguments[0];this.x=e.x,this.y=e.y,this.w=1}else if(2===arguments.length){if("number"==typeof arguments[0]&&"number"==typeof arguments[1]){var n=arguments[0],i=arguments[1];this.x=n,this.y=i,this.w=1}else if(arguments[0]instanceof t&&arguments[1]instanceof t){var r=arguments[0],o=arguments[1];this.x=r.y*o.w-o.y*r.w,this.y=o.x*r.w-r.x*o.w,this.w=r.x*o.y-o.x*r.y}else if(arguments[0]instanceof C&&arguments[1]instanceof C){var s=arguments[0],a=arguments[1];this.x=s.y-a.y,this.y=a.x-s.x,this.w=s.x*a.y-a.x*s.y}}else if(3===arguments.length){var u=arguments[0],l=arguments[1],c=arguments[2];this.x=u,this.y=l,this.w=c}else if(4===arguments.length){var p=arguments[0],h=arguments[1],f=arguments[2],g=arguments[3],d=p.y-h.y,y=h.x-p.x,_=p.x*h.y-h.x*p.y,m=f.y-g.y,v=g.x-f.x,I=f.x*g.y-g.x*f.y;this.x=y*I-v*_,this.y=m*_-d*I,this.w=d*v-m*y}};k.prototype.getY=function(){var t=this.y/this.w;if(v.isNaN(t)||v.isInfinite(t))throw new X;return t},k.prototype.getX=function(){var t=this.x/this.w;if(v.isNaN(t)||v.isInfinite(t))throw new X;return t},k.prototype.getCoordinate=function(){var t=new C;return t.x=this.getX(),t.y=this.getY(),t},k.prototype.interfaces_=function(){return[]},k.prototype.getClass=function(){return k},k.intersection=function(t,e,n,i){var r=t.y-e.y,o=e.x-t.x,s=t.x*e.y-e.x*t.y,a=n.y-i.y,u=i.x-n.x,l=n.x*i.y-i.x*n.y,c=r*u-a*o,p=(o*l-u*s)/c,h=(a*s-r*l)/c;if(v.isNaN(p)||v.isInfinite(p)||v.isNaN(h)||v.isInfinite(h))throw new X;return new C(p,h)};var j=function t(){if(this._minx=null,this._maxx=null,this._miny=null,this._maxy=null,0===arguments.length)this.init();else if(1===arguments.length){if(arguments[0]instanceof C){var e=arguments[0];this.init(e.x,e.x,e.y,e.y)}else if(arguments[0]instanceof t){var n=arguments[0];this.init(n)}}else if(2===arguments.length){var i=arguments[0],r=arguments[1];this.init(i.x,r.x,i.y,r.y)}else if(4===arguments.length){var o=arguments[0],s=arguments[1],a=arguments[2],u=arguments[3];this.init(o,s,a,u)}},H={serialVersionUID:{configurable:!0}};j.prototype.getArea=function(){return this.getWidth()*this.getHeight()},j.prototype.equals=function(t){if(!(t instanceof j))return!1;var e=t;return this.isNull()?e.isNull():this._maxx===e.getMaxX()&&this._maxy===e.getMaxY()&&this._minx===e.getMinX()&&this._miny===e.getMinY()},j.prototype.intersection=function(t){if(this.isNull()||t.isNull()||!this.intersects(t))return new j;var e=this._minx>t._minx?this._minx:t._minx,n=this._miny>t._miny?this._miny:t._miny,i=this._maxx<t._maxx?this._maxx:t._maxx,r=this._maxy<t._maxy?this._maxy:t._maxy;return new j(e,i,n,r)},j.prototype.isNull=function(){return this._maxx<this._minx},j.prototype.getMaxX=function(){return this._maxx},j.prototype.covers=function(){if(1===arguments.length){if(arguments[0]instanceof C){var t=arguments[0];return this.covers(t.x,t.y)}if(arguments[0]instanceof j){var e=arguments[0];return!this.isNull()&&!e.isNull()&&(e.getMinX()>=this._minx&&e.getMaxX()<=this._maxx&&e.getMinY()>=this._miny&&e.getMaxY()<=this._maxy)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];return!this.isNull()&&(n>=this._minx&&n<=this._maxx&&i>=this._miny&&i<=this._maxy)}},j.prototype.intersects=function(){if(1===arguments.length){if(arguments[0]instanceof j){var t=arguments[0];return!this.isNull()&&!t.isNull()&&!(t._minx>this._maxx||t._maxx<this._minx||t._miny>this._maxy||t._maxy<this._miny)}if(arguments[0]instanceof C){var e=arguments[0];return this.intersects(e.x,e.y)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];return!this.isNull()&&!(n>this._maxx||n<this._minx||i>this._maxy||i<this._miny)}},j.prototype.getMinY=function(){return this._miny},j.prototype.getMinX=function(){return this._minx},j.prototype.expandToInclude=function(){if(1===arguments.length){if(arguments[0]instanceof C){var t=arguments[0];this.expandToInclude(t.x,t.y)}else if(arguments[0]instanceof j){var e=arguments[0];if(e.isNull())return null;this.isNull()?(this._minx=e.getMinX(),this._maxx=e.getMaxX(),this._miny=e.getMinY(),this._maxy=e.getMaxY()):(e._minx<this._minx&&(this._minx=e._minx),e._maxx>this._maxx&&(this._maxx=e._maxx),e._miny<this._miny&&(this._miny=e._miny),e._maxy>this._maxy&&(this._maxy=e._maxy))}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this.isNull()?(this._minx=n,this._maxx=n,this._miny=i,this._maxy=i):(n<this._minx&&(this._minx=n),n>this._maxx&&(this._maxx=n),i<this._miny&&(this._miny=i),i>this._maxy&&(this._maxy=i))}},j.prototype.minExtent=function(){if(this.isNull())return 0;var t=this.getWidth(),e=this.getHeight();return t<e?t:e},j.prototype.getWidth=function(){return this.isNull()?0:this._maxx-this._minx},j.prototype.compareTo=function(t){var e=t;return this.isNull()?e.isNull()?0:-1:e.isNull()?1:this._minx<e._minx?-1:this._minx>e._minx?1:this._miny<e._miny?-1:this._miny>e._miny?1:this._maxx<e._maxx?-1:this._maxx>e._maxx?1:this._maxy<e._maxy?-1:this._maxy>e._maxy?1:0},j.prototype.translate=function(t,e){if(this.isNull())return null;this.init(this.getMinX()+t,this.getMaxX()+t,this.getMinY()+e,this.getMaxY()+e)},j.prototype.toString=function(){return"Env["+this._minx+" : "+this._maxx+", "+this._miny+" : "+this._maxy+"]"},j.prototype.setToNull=function(){this._minx=0,this._maxx=-1,this._miny=0,this._maxy=-1},j.prototype.getHeight=function(){return this.isNull()?0:this._maxy-this._miny},j.prototype.maxExtent=function(){if(this.isNull())return 0;var t=this.getWidth(),e=this.getHeight();return t>e?t:e},j.prototype.expandBy=function(){if(1===arguments.length){var t=arguments[0];this.expandBy(t,t)}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(this.isNull())return null;this._minx-=e,this._maxx+=e,this._miny-=n,this._maxy+=n,(this._minx>this._maxx||this._miny>this._maxy)&&this.setToNull()}},j.prototype.contains=function(){if(1===arguments.length){if(arguments[0]instanceof j){var t=arguments[0];return this.covers(t)}if(arguments[0]instanceof C){var e=arguments[0];return this.covers(e)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];return this.covers(n,i)}},j.prototype.centre=function(){return this.isNull()?null:new C((this.getMinX()+this.getMaxX())/2,(this.getMinY()+this.getMaxY())/2)},j.prototype.init=function(){if(0===arguments.length)this.setToNull();else if(1===arguments.length){if(arguments[0]instanceof C){var t=arguments[0];this.init(t.x,t.x,t.y,t.y)}else if(arguments[0]instanceof j){var e=arguments[0];this._minx=e._minx,this._maxx=e._maxx,this._miny=e._miny,this._maxy=e._maxy}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this.init(n.x,i.x,n.y,i.y)}else if(4===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2],a=arguments[3];r<o?(this._minx=r,this._maxx=o):(this._minx=o,this._maxx=r),s<a?(this._miny=s,this._maxy=a):(this._miny=a,this._maxy=s)}},j.prototype.getMaxY=function(){return this._maxy},j.prototype.distance=function(t){if(this.intersects(t))return 0;var e=0;this._maxx<t._minx?e=t._minx-this._maxx:this._minx>t._maxx&&(e=this._minx-t._maxx);var n=0;return this._maxy<t._miny?n=t._miny-this._maxy:this._miny>t._maxy&&(n=this._miny-t._maxy),0===e?n:0===n?e:Math.sqrt(e*e+n*n)},j.prototype.hashCode=function(){var t=17;return t=37*t+C.hashCode(this._minx),t=37*t+C.hashCode(this._maxx),t=37*t+C.hashCode(this._miny),t=37*t+C.hashCode(this._maxy)},j.prototype.interfaces_=function(){return[E,e]},j.prototype.getClass=function(){return j},j.intersects=function(){if(3===arguments.length){var t=arguments[0],e=arguments[1],n=arguments[2];return n.x>=(t.x<e.x?t.x:e.x)&&n.x<=(t.x>e.x?t.x:e.x)&&n.y>=(t.y<e.y?t.y:e.y)&&n.y<=(t.y>e.y?t.y:e.y)}if(4===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2],s=arguments[3],a=Math.min(o.x,s.x),u=Math.max(o.x,s.x),l=Math.min(i.x,r.x),c=Math.max(i.x,r.x);return!(l>u)&&(!(c<a)&&(a=Math.min(o.y,s.y),u=Math.max(o.y,s.y),l=Math.min(i.y,r.y),c=Math.max(i.y,r.y),!(l>u)&&!(c<a)))}},H.serialVersionUID.get=function(){return 0x51845cd552189800},Object.defineProperties(j,H);var W={typeStr:/^\s*(\w+)\s*\(\s*(.*)\s*\)\s*$/,emptyTypeStr:/^\s*(\w+)\s*EMPTY\s*$/,spaces:/\s+/,parenComma:/\)\s*,\s*\(/,doubleParenComma:/\)\s*\)\s*,\s*\(\s*\(/,trimParens:/^\s*\(?(.*?)\)?\s*$/},K=function(t){this.geometryFactory=t||new _e};K.prototype.read=function(t){var e,n,i;t=t.replace(/[\n\r]/g," ");var r=W.typeStr.exec(t);if(-1!==t.search("EMPTY")&&((r=W.emptyTypeStr.exec(t))[2]=void 0),r&&(n=r[1].toLowerCase(),i=r[2],Q[n]&&(e=Q[n].apply(this,[i]))),void 0===e)throw new Error("Could not parse WKT "+t);return e},K.prototype.write=function(t){return this.extractGeometry(t)},K.prototype.extractGeometry=function(t){var e=t.getGeometryType().toLowerCase();if(!J[e])return null;var n=e.toUpperCase();return t.isEmpty()?n+" EMPTY":n+"("+J[e].apply(this,[t])+")"};var J={coordinate:function(t){return t.x+" "+t.y},point:function(t){return J.coordinate.call(this,t._coordinates._coordinates[0])},multipoint:function(t){for(var e=[],n=0,i=t._geometries.length;n<i;++n)e.push("("+J.point.apply(this,[t._geometries[n]])+")");return e.join(",")},linestring:function(t){for(var e=[],n=0,i=t._points._coordinates.length;n<i;++n)e.push(J.coordinate.apply(this,[t._points._coordinates[n]]));return e.join(",")},linearring:function(t){for(var e=[],n=0,i=t._points._coordinates.length;n<i;++n)e.push(J.coordinate.apply(this,[t._points._coordinates[n]]));return e.join(",")},multilinestring:function(t){for(var e=[],n=0,i=t._geometries.length;n<i;++n)e.push("("+J.linestring.apply(this,[t._geometries[n]])+")");return e.join(",")},polygon:function(t){var e=[];e.push("("+J.linestring.apply(this,[t._shell])+")");for(var n=0,i=t._holes.length;n<i;++n)e.push("("+J.linestring.apply(this,[t._holes[n]])+")");return e.join(",")},multipolygon:function(t){for(var e=[],n=0,i=t._geometries.length;n<i;++n)e.push("("+J.polygon.apply(this,[t._geometries[n]])+")");return e.join(",")},geometrycollection:function(t){for(var e=[],n=0,i=t._geometries.length;n<i;++n)e.push(this.extractGeometry(t._geometries[n]));return e.join(",")}},Q={point:function(t){if(void 0===t)return this.geometryFactory.createPoint();var e=t.trim().split(W.spaces);return this.geometryFactory.createPoint(new C(Number.parseFloat(e[0]),Number.parseFloat(e[1])))},multipoint:function(t){if(void 0===t)return this.geometryFactory.createMultiPoint();for(var e,n=t.trim().split(","),i=[],r=0,o=n.length;r<o;++r)e=n[r].replace(W.trimParens,"$1"),i.push(Q.point.apply(this,[e]));return this.geometryFactory.createMultiPoint(i)},linestring:function(t){if(void 0===t)return this.geometryFactory.createLineString();for(var e,n=t.trim().split(","),i=[],r=0,o=n.length;r<o;++r)e=n[r].trim().split(W.spaces),i.push(new C(Number.parseFloat(e[0]),Number.parseFloat(e[1])));return this.geometryFactory.createLineString(i)},linearring:function(t){if(void 0===t)return this.geometryFactory.createLinearRing();for(var e,n=t.trim().split(","),i=[],r=0,o=n.length;r<o;++r)e=n[r].trim().split(W.spaces),i.push(new C(Number.parseFloat(e[0]),Number.parseFloat(e[1])));return this.geometryFactory.createLinearRing(i)},multilinestring:function(t){if(void 0===t)return this.geometryFactory.createMultiLineString();for(var e,n=t.trim().split(W.parenComma),i=[],r=0,o=n.length;r<o;++r)e=n[r].replace(W.trimParens,"$1"),i.push(Q.linestring.apply(this,[e]));return this.geometryFactory.createMultiLineString(i)},polygon:function(t){if(void 0===t)return this.geometryFactory.createPolygon();for(var e,n,i,r,o=t.trim().split(W.parenComma),s=[],a=0,u=o.length;a<u;++a)e=o[a].replace(W.trimParens,"$1"),n=Q.linestring.apply(this,[e]),i=this.geometryFactory.createLinearRing(n._points),0===a?r=i:s.push(i);return this.geometryFactory.createPolygon(r,s)},multipolygon:function(t){if(void 0===t)return this.geometryFactory.createMultiPolygon();for(var e,n=t.trim().split(W.doubleParenComma),i=[],r=0,o=n.length;r<o;++r)e=n[r].replace(W.trimParens,"$1"),i.push(Q.polygon.apply(this,[e]));return this.geometryFactory.createMultiPolygon(i)},geometrycollection:function(t){if(void 0===t)return this.geometryFactory.createGeometryCollection();for(var e=(t=t.replace(/,\s*([A-Za-z])/g,"|$1")).trim().split("|"),n=[],i=0,r=e.length;i<r;++i)n.push(this.read(e[i]));return this.geometryFactory.createGeometryCollection(n)}},Z=function(t){this.parser=new K(t)};Z.prototype.write=function(t){return this.parser.write(t)},Z.toLineString=function(t,e){if(2!==arguments.length)throw new Error("Not implemented");return"LINESTRING ( "+t.x+" "+t.y+", "+e.x+" "+e.y+" )"};var $=function(t){function e(e){t.call(this,e),this.name="RuntimeException",this.message=e,this.stack=(new t).stack}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e}(Error),tt=function(t){function e(){if(t.call(this),0===arguments.length)t.call(this);else if(1===arguments.length){var e=arguments[0];t.call(this,e)}}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}($),et=function(){};et.prototype.interfaces_=function(){return[]},et.prototype.getClass=function(){return et},et.shouldNeverReachHere=function(){if(0===arguments.length)et.shouldNeverReachHere(null);else if(1===arguments.length){var t=arguments[0];throw new tt("Should never reach here"+(null!==t?": "+t:""))}},et.isTrue=function(){var t,e;if(1===arguments.length)t=arguments[0],et.isTrue(t,null);else if(2===arguments.length&&(t=arguments[0],e=arguments[1],!t))throw null===e?new tt:new tt(e)},et.equals=function(){var t,e,n;if(2===arguments.length)t=arguments[0],e=arguments[1],et.equals(t,e,null);else if(3===arguments.length&&(t=arguments[0],e=arguments[1],n=arguments[2],!e.equals(t)))throw new tt("Expected "+t+" but encountered "+e+(null!==n?": "+n:""))};var nt=function(){this._result=null,this._inputLines=Array(2).fill().map(function(){return Array(2)}),this._intPt=new Array(2).fill(null),this._intLineIndex=null,this._isProper=null,this._pa=null,this._pb=null,this._precisionModel=null,this._intPt[0]=new C,this._intPt[1]=new C,this._pa=this._intPt[0],this._pb=this._intPt[1],this._result=0},it={DONT_INTERSECT:{configurable:!0},DO_INTERSECT:{configurable:!0},COLLINEAR:{configurable:!0},NO_INTERSECTION:{configurable:!0},POINT_INTERSECTION:{configurable:!0},COLLINEAR_INTERSECTION:{configurable:!0}};nt.prototype.getIndexAlongSegment=function(t,e){return this.computeIntLineIndex(),this._intLineIndex[t][e]},nt.prototype.getTopologySummary=function(){var t=new D;return this.isEndPoint()&&t.append(" endpoint"),this._isProper&&t.append(" proper"),this.isCollinear()&&t.append(" collinear"),t.toString()},nt.prototype.computeIntersection=function(t,e,n,i){this._inputLines[0][0]=t,this._inputLines[0][1]=e,this._inputLines[1][0]=n,this._inputLines[1][1]=i,this._result=this.computeIntersect(t,e,n,i)},nt.prototype.getIntersectionNum=function(){return this._result},nt.prototype.computeIntLineIndex=function(){if(0===arguments.length)null===this._intLineIndex&&(this._intLineIndex=Array(2).fill().map(function(){return Array(2)}),this.computeIntLineIndex(0),this.computeIntLineIndex(1));else if(1===arguments.length){var t=arguments[0];this.getEdgeDistance(t,0)>this.getEdgeDistance(t,1)?(this._intLineIndex[t][0]=0,this._intLineIndex[t][1]=1):(this._intLineIndex[t][0]=1,this._intLineIndex[t][1]=0)}},nt.prototype.isProper=function(){return this.hasIntersection()&&this._isProper},nt.prototype.setPrecisionModel=function(t){this._precisionModel=t},nt.prototype.isInteriorIntersection=function(){if(0===arguments.length)return!!this.isInteriorIntersection(0)||!!this.isInteriorIntersection(1);if(1===arguments.length){for(var t=arguments[0],e=0;e<this._result;e++)if(!this._intPt[e].equals2D(this._inputLines[t][0])&&!this._intPt[e].equals2D(this._inputLines[t][1]))return!0;return!1}},nt.prototype.getIntersection=function(t){return this._intPt[t]},nt.prototype.isEndPoint=function(){return this.hasIntersection()&&!this._isProper},nt.prototype.hasIntersection=function(){return this._result!==nt.NO_INTERSECTION},nt.prototype.getEdgeDistance=function(t,e){return nt.computeEdgeDistance(this._intPt[e],this._inputLines[t][0],this._inputLines[t][1])},nt.prototype.isCollinear=function(){return this._result===nt.COLLINEAR_INTERSECTION},nt.prototype.toString=function(){return Z.toLineString(this._inputLines[0][0],this._inputLines[0][1])+" - "+Z.toLineString(this._inputLines[1][0],this._inputLines[1][1])+this.getTopologySummary()},nt.prototype.getEndpoint=function(t,e){return this._inputLines[t][e]},nt.prototype.isIntersection=function(t){for(var e=0;e<this._result;e++)if(this._intPt[e].equals2D(t))return!0;return!1},nt.prototype.getIntersectionAlongSegment=function(t,e){return this.computeIntLineIndex(),this._intPt[this._intLineIndex[t][e]]},nt.prototype.interfaces_=function(){return[]},nt.prototype.getClass=function(){return nt},nt.computeEdgeDistance=function(t,e,n){var i=Math.abs(n.x-e.x),r=Math.abs(n.y-e.y),o=-1;if(t.equals(e))o=0;else if(t.equals(n))o=i>r?i:r;else{var s=Math.abs(t.x-e.x),a=Math.abs(t.y-e.y);0!==(o=i>r?s:a)||t.equals(e)||(o=Math.max(s,a))}return et.isTrue(!(0===o&&!t.equals(e)),"Bad distance calculation"),o},nt.nonRobustComputeEdgeDistance=function(t,e,n){var i=t.x-e.x,r=t.y-e.y,o=Math.sqrt(i*i+r*r);return et.isTrue(!(0===o&&!t.equals(e)),"Invalid distance calculation"),o},it.DONT_INTERSECT.get=function(){return 0},it.DO_INTERSECT.get=function(){return 1},it.COLLINEAR.get=function(){return 2},it.NO_INTERSECTION.get=function(){return 0},it.POINT_INTERSECTION.get=function(){return 1},it.COLLINEAR_INTERSECTION.get=function(){return 2},Object.defineProperties(nt,it);var rt=function(t){function e(){t.apply(this,arguments)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.isInSegmentEnvelopes=function(t){var e=new j(this._inputLines[0][0],this._inputLines[0][1]),n=new j(this._inputLines[1][0],this._inputLines[1][1]);return e.contains(t)&&n.contains(t)},e.prototype.computeIntersection=function(){if(3!==arguments.length)return t.prototype.computeIntersection.apply(this,arguments);var e=arguments[0],n=arguments[1],i=arguments[2];if(this._isProper=!1,j.intersects(n,i,e)&&0===at.orientationIndex(n,i,e)&&0===at.orientationIndex(i,n,e))return this._isProper=!0,(e.equals(n)||e.equals(i))&&(this._isProper=!1),this._result=t.POINT_INTERSECTION,null;this._result=t.NO_INTERSECTION},e.prototype.normalizeToMinimum=function(t,e,n,i,r){r.x=this.smallestInAbsValue(t.x,e.x,n.x,i.x),r.y=this.smallestInAbsValue(t.y,e.y,n.y,i.y),t.x-=r.x,t.y-=r.y,e.x-=r.x,e.y-=r.y,n.x-=r.x,n.y-=r.y,i.x-=r.x,i.y-=r.y},e.prototype.safeHCoordinateIntersection=function(t,n,i,r){var o=null;try{o=k.intersection(t,n,i,r)}catch(s){if(!(s instanceof X))throw s;o=e.nearestEndpoint(t,n,i,r)}return o},e.prototype.intersection=function(t,n,i,r){var o=this.intersectionWithNormalization(t,n,i,r);return this.isInSegmentEnvelopes(o)||(o=new C(e.nearestEndpoint(t,n,i,r))),null!==this._precisionModel&&this._precisionModel.makePrecise(o),o},e.prototype.smallestInAbsValue=function(t,e,n,i){var r=t,o=Math.abs(r);return Math.abs(e)<o&&(r=e,o=Math.abs(e)),Math.abs(n)<o&&(r=n,o=Math.abs(n)),Math.abs(i)<o&&(r=i),r},e.prototype.checkDD=function(t,e,n,i,r){var o=q.intersection(t,e,n,i),s=this.isInSegmentEnvelopes(o);Y.out.println("DD in env = "+s+"  --------------------- "+o),r.distance(o)>1e-4&&Y.out.println("Distance = "+r.distance(o))},e.prototype.intersectionWithNormalization=function(t,e,n,i){var r=new C(t),o=new C(e),s=new C(n),a=new C(i),u=new C;this.normalizeToEnvCentre(r,o,s,a,u);var l=this.safeHCoordinateIntersection(r,o,s,a);return l.x+=u.x,l.y+=u.y,l},e.prototype.computeCollinearIntersection=function(e,n,i,r){var o=j.intersects(e,n,i),s=j.intersects(e,n,r),a=j.intersects(i,r,e),u=j.intersects(i,r,n);return o&&s?(this._intPt[0]=i,this._intPt[1]=r,t.COLLINEAR_INTERSECTION):a&&u?(this._intPt[0]=e,this._intPt[1]=n,t.COLLINEAR_INTERSECTION):o&&a?(this._intPt[0]=i,this._intPt[1]=e,!i.equals(e)||s||u?t.COLLINEAR_INTERSECTION:t.POINT_INTERSECTION):o&&u?(this._intPt[0]=i,this._intPt[1]=n,!i.equals(n)||s||a?t.COLLINEAR_INTERSECTION:t.POINT_INTERSECTION):s&&a?(this._intPt[0]=r,this._intPt[1]=e,!r.equals(e)||o||u?t.COLLINEAR_INTERSECTION:t.POINT_INTERSECTION):s&&u?(this._intPt[0]=r,this._intPt[1]=n,!r.equals(n)||o||a?t.COLLINEAR_INTERSECTION:t.POINT_INTERSECTION):t.NO_INTERSECTION},e.prototype.normalizeToEnvCentre=function(t,e,n,i,r){var o=t.x<e.x?t.x:e.x,s=t.y<e.y?t.y:e.y,a=t.x>e.x?t.x:e.x,u=t.y>e.y?t.y:e.y,l=n.x<i.x?n.x:i.x,c=n.y<i.y?n.y:i.y,p=n.x>i.x?n.x:i.x,h=n.y>i.y?n.y:i.y,f=((o>l?o:l)+(a<p?a:p))/2,g=((s>c?s:c)+(u<h?u:h))/2;r.x=f,r.y=g,t.x-=r.x,t.y-=r.y,e.x-=r.x,e.y-=r.y,n.x-=r.x,n.y-=r.y,i.x-=r.x,i.y-=r.y},e.prototype.computeIntersect=function(e,n,i,r){if(this._isProper=!1,!j.intersects(e,n,i,r))return t.NO_INTERSECTION;var o=at.orientationIndex(e,n,i),s=at.orientationIndex(e,n,r);if(o>0&&s>0||o<0&&s<0)return t.NO_INTERSECTION;var a=at.orientationIndex(i,r,e),u=at.orientationIndex(i,r,n);if(a>0&&u>0||a<0&&u<0)return t.NO_INTERSECTION;return 0===o&&0===s&&0===a&&0===u?this.computeCollinearIntersection(e,n,i,r):(0===o||0===s||0===a||0===u?(this._isProper=!1,e.equals2D(i)||e.equals2D(r)?this._intPt[0]=e:n.equals2D(i)||n.equals2D(r)?this._intPt[0]=n:0===o?this._intPt[0]=new C(i):0===s?this._intPt[0]=new C(r):0===a?this._intPt[0]=new C(e):0===u&&(this._intPt[0]=new C(n))):(this._isProper=!0,this._intPt[0]=this.intersection(e,n,i,r)),t.POINT_INTERSECTION)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e.nearestEndpoint=function(t,e,n,i){var r=t,o=at.distancePointLine(t,n,i),s=at.distancePointLine(e,n,i);return s<o&&(o=s,r=e),(s=at.distancePointLine(n,t,e))<o&&(o=s,r=n),(s=at.distancePointLine(i,t,e))<o&&(o=s,r=i),r},e}(nt),ot=function(){};ot.prototype.interfaces_=function(){return[]},ot.prototype.getClass=function(){return ot},ot.orientationIndex=function(t,e,n){var i=e.x-t.x,r=e.y-t.y,o=n.x-e.x,s=n.y-e.y;return ot.signOfDet2x2(i,r,o,s)},ot.signOfDet2x2=function(t,e,n,i){var r=null,o=null,s=null;if(r=1,0===t||0===i)return 0===e||0===n?0:e>0?n>0?-r:r:n>0?r:-r;if(0===e||0===n)return i>0?t>0?r:-r:t>0?-r:r;if(e>0?i>0?e<=i||(r=-r,o=t,t=n,n=o,o=e,e=i,i=o):e<=-i?(r=-r,n=-n,i=-i):(o=t,t=-n,n=o,o=e,e=-i,i=o):i>0?-e<=i?(r=-r,t=-t,e=-e):(o=-t,t=n,n=o,o=-e,e=i,i=o):e>=i?(t=-t,e=-e,n=-n,i=-i):(r=-r,o=-t,t=-n,n=o,o=-e,e=-i,i=o),t>0){if(!(n>0))return r;if(!(t<=n))return r}else{if(n>0)return-r;if(!(t>=n))return-r;r=-r,t=-t,n=-n}for(;;){if(s=Math.floor(n/t),n-=s*t,(i-=s*e)<0)return-r;if(i>e)return r;if(t>n+n){if(e<i+i)return r}else{if(e>i+i)return-r;n=t-n,i=e-i,r=-r}if(0===i)return 0===n?0:-r;if(0===n)return r;if(s=Math.floor(t/n),t-=s*n,(e-=s*i)<0)return r;if(e>i)return-r;if(n>t+t){if(i<e+e)return-r}else{if(i>e+e)return r;t=n-t,e=i-e,r=-r}if(0===e)return 0===t?0:r;if(0===t)return-r}};var st=function(){this._p=null,this._crossingCount=0,this._isPointOnSegment=!1;var t=arguments[0];this._p=t};st.prototype.countSegment=function(t,e){if(t.x<this._p.x&&e.x<this._p.x)return null;if(this._p.x===e.x&&this._p.y===e.y)return this._isPointOnSegment=!0,null;if(t.y===this._p.y&&e.y===this._p.y){var n=t.x,i=e.x;return n>i&&(n=e.x,i=t.x),this._p.x>=n&&this._p.x<=i&&(this._isPointOnSegment=!0),null}if(t.y>this._p.y&&e.y<=this._p.y||e.y>this._p.y&&t.y<=this._p.y){var r=t.x-this._p.x,o=t.y-this._p.y,s=e.x-this._p.x,a=e.y-this._p.y,u=ot.signOfDet2x2(r,o,s,a);if(0===u)return this._isPointOnSegment=!0,null;a<o&&(u=-u),u>0&&this._crossingCount++}},st.prototype.isPointInPolygon=function(){return this.getLocation()!==w.EXTERIOR},st.prototype.getLocation=function(){return this._isPointOnSegment?w.BOUNDARY:this._crossingCount%2==1?w.INTERIOR:w.EXTERIOR},st.prototype.isOnSegment=function(){return this._isPointOnSegment},st.prototype.interfaces_=function(){return[]},st.prototype.getClass=function(){return st},st.locatePointInRing=function(){if(arguments[0]instanceof C&&T(arguments[1],V)){for(var t=arguments[0],e=arguments[1],n=new st(t),i=new C,r=new C,o=1;o<e.size();o++)if(e.getCoordinate(o,i),e.getCoordinate(o-1,r),n.countSegment(i,r),n.isOnSegment())return n.getLocation();return n.getLocation()}if(arguments[0]instanceof C&&arguments[1]instanceof Array){for(var s=arguments[0],a=arguments[1],u=new st(s),l=1;l<a.length;l++){var c=a[l],p=a[l-1];if(u.countSegment(c,p),u.isOnSegment())return u.getLocation()}return u.getLocation()}};var at=function(){},ut={CLOCKWISE:{configurable:!0},RIGHT:{configurable:!0},COUNTERCLOCKWISE:{configurable:!0},LEFT:{configurable:!0},COLLINEAR:{configurable:!0},STRAIGHT:{configurable:!0}};at.prototype.interfaces_=function(){return[]},at.prototype.getClass=function(){return at},at.orientationIndex=function(t,e,n){return q.orientationIndex(t,e,n)},at.signedArea=function(){if(arguments[0]instanceof Array){var t=arguments[0];if(t.length<3)return 0;for(var e=0,n=t[0].x,i=1;i<t.length-1;i++){var r=t[i].x-n,o=t[i+1].y;e+=r*(t[i-1].y-o)}return e/2}if(T(arguments[0],V)){var s=arguments[0],a=s.size();if(a<3)return 0;var u=new C,l=new C,c=new C;s.getCoordinate(0,l),s.getCoordinate(1,c);var p=l.x;c.x-=p;for(var h=0,f=1;f<a-1;f++)u.y=l.y,l.x=c.x,l.y=c.y,s.getCoordinate(f+1,c),c.x-=p,h+=l.x*(u.y-c.y);return h/2}},at.distanceLineLine=function(t,e,n,i){if(t.equals(e))return at.distancePointLine(t,n,i);if(n.equals(i))return at.distancePointLine(i,t,e);var r=!1;if(j.intersects(t,e,n,i)){var o=(e.x-t.x)*(i.y-n.y)-(e.y-t.y)*(i.x-n.x);if(0===o)r=!0;else{var s=(t.y-n.y)*(i.x-n.x)-(t.x-n.x)*(i.y-n.y),a=((t.y-n.y)*(e.x-t.x)-(t.x-n.x)*(e.y-t.y))/o,u=s/o;(u<0||u>1||a<0||a>1)&&(r=!0)}}else r=!0;return r?R.min(at.distancePointLine(t,n,i),at.distancePointLine(e,n,i),at.distancePointLine(n,t,e),at.distancePointLine(i,t,e)):0},at.isPointInRing=function(t,e){return at.locatePointInRing(t,e)!==w.EXTERIOR},at.computeLength=function(t){var e=t.size();if(e<=1)return 0;var n=0,i=new C;t.getCoordinate(0,i);for(var r=i.x,o=i.y,s=1;s<e;s++){t.getCoordinate(s,i);var a=i.x,u=i.y,l=a-r,c=u-o;n+=Math.sqrt(l*l+c*c),r=a,o=u}return n},at.isCCW=function(t){var e=t.length-1;if(e<3)throw new m("Ring has fewer than 4 points, so orientation cannot be determined");for(var n=t[0],i=0,r=1;r<=e;r++){var o=t[r];o.y>n.y&&(n=o,i=r)}var s=i;do{(s-=1)<0&&(s=e)}while(t[s].equals2D(n)&&s!==i);var a=i;do{a=(a+1)%e}while(t[a].equals2D(n)&&a!==i);var u=t[s],l=t[a];if(u.equals2D(n)||l.equals2D(n)||u.equals2D(l))return!1;var c=at.computeOrientation(u,n,l),p=!1;return p=0===c?u.x>l.x:c>0,p},at.locatePointInRing=function(t,e){return st.locatePointInRing(t,e)},at.distancePointLinePerpendicular=function(t,e,n){var i=(n.x-e.x)*(n.x-e.x)+(n.y-e.y)*(n.y-e.y),r=((e.y-t.y)*(n.x-e.x)-(e.x-t.x)*(n.y-e.y))/i;return Math.abs(r)*Math.sqrt(i)},at.computeOrientation=function(t,e,n){return at.orientationIndex(t,e,n)},at.distancePointLine=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];if(0===e.length)throw new m("Line array must contain at least one vertex");for(var n=t.distance(e[0]),i=0;i<e.length-1;i++){var r=at.distancePointLine(t,e[i],e[i+1]);r<n&&(n=r)}return n}if(3===arguments.length){var o=arguments[0],s=arguments[1],a=arguments[2];if(s.x===a.x&&s.y===a.y)return o.distance(s);var u=(a.x-s.x)*(a.x-s.x)+(a.y-s.y)*(a.y-s.y),l=((o.x-s.x)*(a.x-s.x)+(o.y-s.y)*(a.y-s.y))/u;if(l<=0)return o.distance(s);if(l>=1)return o.distance(a);var c=((s.y-o.y)*(a.x-s.x)-(s.x-o.x)*(a.y-s.y))/u;return Math.abs(c)*Math.sqrt(u)}},at.isOnLine=function(t,e){for(var n=new rt,i=1;i<e.length;i++){var r=e[i-1],o=e[i];if(n.computeIntersection(t,r,o),n.hasIntersection())return!0}return!1},ut.CLOCKWISE.get=function(){return-1},ut.RIGHT.get=function(){return at.CLOCKWISE},ut.COUNTERCLOCKWISE.get=function(){return 1},ut.LEFT.get=function(){return at.COUNTERCLOCKWISE},ut.COLLINEAR.get=function(){return 0},ut.STRAIGHT.get=function(){return at.COLLINEAR},Object.defineProperties(at,ut);var lt=function(){};lt.prototype.filter=function(t){},lt.prototype.interfaces_=function(){return[]},lt.prototype.getClass=function(){return lt};var ct=function(){var t=arguments[0];this._envelope=null,this._factory=null,this._SRID=null,this._userData=null,this._factory=t,this._SRID=t.getSRID()},pt={serialVersionUID:{configurable:!0},SORTINDEX_POINT:{configurable:!0},SORTINDEX_MULTIPOINT:{configurable:!0},SORTINDEX_LINESTRING:{configurable:!0},SORTINDEX_LINEARRING:{configurable:!0},SORTINDEX_MULTILINESTRING:{configurable:!0},SORTINDEX_POLYGON:{configurable:!0},SORTINDEX_MULTIPOLYGON:{configurable:!0},SORTINDEX_GEOMETRYCOLLECTION:{configurable:!0},geometryChangedFilter:{configurable:!0}};ct.prototype.isGeometryCollection=function(){return this.getSortIndex()===ct.SORTINDEX_GEOMETRYCOLLECTION},ct.prototype.getFactory=function(){return this._factory},ct.prototype.getGeometryN=function(t){return this},ct.prototype.getArea=function(){return 0},ct.prototype.isRectangle=function(){return!1},ct.prototype.equals=function(){if(arguments[0]instanceof ct){var t=arguments[0];return null!==t&&this.equalsTopo(t)}if(arguments[0]instanceof Object){var e=arguments[0];if(!(e instanceof ct))return!1;var n=e;return this.equalsExact(n)}},ct.prototype.equalsExact=function(t){return this===t||this.equalsExact(t,0)},ct.prototype.geometryChanged=function(){this.apply(ct.geometryChangedFilter)},ct.prototype.geometryChangedAction=function(){this._envelope=null},ct.prototype.equalsNorm=function(t){return null!==t&&this.norm().equalsExact(t.norm())},ct.prototype.getLength=function(){return 0},ct.prototype.getNumGeometries=function(){return 1},ct.prototype.compareTo=function(){if(1===arguments.length){var t=arguments[0],e=t;return this.getSortIndex()!==e.getSortIndex()?this.getSortIndex()-e.getSortIndex():this.isEmpty()&&e.isEmpty()?0:this.isEmpty()?-1:e.isEmpty()?1:this.compareToSameClass(t)}if(2===arguments.length){var n=arguments[0],i=arguments[1];return this.getSortIndex()!==n.getSortIndex()?this.getSortIndex()-n.getSortIndex():this.isEmpty()&&n.isEmpty()?0:this.isEmpty()?-1:n.isEmpty()?1:this.compareToSameClass(n,i)}},ct.prototype.getUserData=function(){return this._userData},ct.prototype.getSRID=function(){return this._SRID},ct.prototype.getEnvelope=function(){return this.getFactory().toGeometry(this.getEnvelopeInternal())},ct.prototype.checkNotGeometryCollection=function(t){if(t.getSortIndex()===ct.SORTINDEX_GEOMETRYCOLLECTION)throw new m("This method does not support GeometryCollection arguments")},ct.prototype.equal=function(t,e,n){return 0===n?t.equals(e):t.distance(e)<=n},ct.prototype.norm=function(){var t=this.copy();return t.normalize(),t},ct.prototype.getPrecisionModel=function(){return this._factory.getPrecisionModel()},ct.prototype.getEnvelopeInternal=function(){return null===this._envelope&&(this._envelope=this.computeEnvelopeInternal()),new j(this._envelope)},ct.prototype.setSRID=function(t){this._SRID=t},ct.prototype.setUserData=function(t){this._userData=t},ct.prototype.compare=function(t,e){for(var n=t.iterator(),i=e.iterator();n.hasNext()&&i.hasNext();){var r=n.next(),o=i.next(),s=r.compareTo(o);if(0!==s)return s}return n.hasNext()?1:i.hasNext()?-1:0},ct.prototype.hashCode=function(){return this.getEnvelopeInternal().hashCode()},ct.prototype.isGeometryCollectionOrDerived=function(){return this.getSortIndex()===ct.SORTINDEX_GEOMETRYCOLLECTION||this.getSortIndex()===ct.SORTINDEX_MULTIPOINT||this.getSortIndex()===ct.SORTINDEX_MULTILINESTRING||this.getSortIndex()===ct.SORTINDEX_MULTIPOLYGON},ct.prototype.interfaces_=function(){return[x,E,e]},ct.prototype.getClass=function(){return ct},ct.hasNonEmptyElements=function(t){for(var e=0;e<t.length;e++)if(!t[e].isEmpty())return!0;return!1},ct.hasNullElements=function(t){for(var e=0;e<t.length;e++)if(null===t[e])return!0;return!1},pt.serialVersionUID.get=function(){return 0x799ea46522854c00},pt.SORTINDEX_POINT.get=function(){return 0},pt.SORTINDEX_MULTIPOINT.get=function(){return 1},pt.SORTINDEX_LINESTRING.get=function(){return 2},pt.SORTINDEX_LINEARRING.get=function(){return 3},pt.SORTINDEX_MULTILINESTRING.get=function(){return 4},pt.SORTINDEX_POLYGON.get=function(){return 5},pt.SORTINDEX_MULTIPOLYGON.get=function(){return 6},pt.SORTINDEX_GEOMETRYCOLLECTION.get=function(){return 7},pt.geometryChangedFilter.get=function(){return ht},Object.defineProperties(ct,pt);var ht=function(){};ht.interfaces_=function(){return[lt]},ht.filter=function(t){t.geometryChangedAction()};var ft=function(){};ft.prototype.filter=function(t){},ft.prototype.interfaces_=function(){return[]},ft.prototype.getClass=function(){return ft};var gt=function(){},dt={Mod2BoundaryNodeRule:{configurable:!0},EndPointBoundaryNodeRule:{configurable:!0},MultiValentEndPointBoundaryNodeRule:{configurable:!0},MonoValentEndPointBoundaryNodeRule:{configurable:!0},MOD2_BOUNDARY_RULE:{configurable:!0},ENDPOINT_BOUNDARY_RULE:{configurable:!0},MULTIVALENT_ENDPOINT_BOUNDARY_RULE:{configurable:!0},MONOVALENT_ENDPOINT_BOUNDARY_RULE:{configurable:!0},OGC_SFS_BOUNDARY_RULE:{configurable:!0}};gt.prototype.isInBoundary=function(t){},gt.prototype.interfaces_=function(){return[]},gt.prototype.getClass=function(){return gt},dt.Mod2BoundaryNodeRule.get=function(){return yt},dt.EndPointBoundaryNodeRule.get=function(){return _t},dt.MultiValentEndPointBoundaryNodeRule.get=function(){return mt},dt.MonoValentEndPointBoundaryNodeRule.get=function(){return vt},dt.MOD2_BOUNDARY_RULE.get=function(){return new yt},dt.ENDPOINT_BOUNDARY_RULE.get=function(){return new _t},dt.MULTIVALENT_ENDPOINT_BOUNDARY_RULE.get=function(){return new mt},dt.MONOVALENT_ENDPOINT_BOUNDARY_RULE.get=function(){return new vt},dt.OGC_SFS_BOUNDARY_RULE.get=function(){return gt.MOD2_BOUNDARY_RULE},Object.defineProperties(gt,dt);var yt=function(){};yt.prototype.isInBoundary=function(t){return t%2==1},yt.prototype.interfaces_=function(){return[gt]},yt.prototype.getClass=function(){return yt};var _t=function(){};_t.prototype.isInBoundary=function(t){return t>0},_t.prototype.interfaces_=function(){return[gt]},_t.prototype.getClass=function(){return _t};var mt=function(){};mt.prototype.isInBoundary=function(t){return t>1},mt.prototype.interfaces_=function(){return[gt]},mt.prototype.getClass=function(){return mt};var vt=function(){};vt.prototype.isInBoundary=function(t){return 1===t},vt.prototype.interfaces_=function(){return[gt]},vt.prototype.getClass=function(){return vt};var It=function(){};It.prototype.add=function(){},It.prototype.addAll=function(){},It.prototype.isEmpty=function(){},It.prototype.iterator=function(){},It.prototype.size=function(){},It.prototype.toArray=function(){},It.prototype.remove=function(){},(n.prototype=new Error).name="IndexOutOfBoundsException";var Et=function(){};Et.prototype.hasNext=function(){},Et.prototype.next=function(){},Et.prototype.remove=function(){};var xt=function(t){function e(){t.apply(this,arguments)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.get=function(){},e.prototype.set=function(){},e.prototype.isEmpty=function(){},e}(It);(i.prototype=new Error).name="NoSuchElementException";var Nt=function(t){function e(){t.call(this),this.array_=[],arguments[0]instanceof It&&this.addAll(arguments[0])}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.ensureCapacity=function(){},e.prototype.interfaces_=function(){return[t,It]},e.prototype.add=function(t){return 1===arguments.length?this.array_.push(t):this.array_.splice(arguments[0],arguments[1]),!0},e.prototype.clear=function(){this.array_=[]},e.prototype.addAll=function(t){for(var e=t.iterator();e.hasNext();)this.add(e.next());return!0},e.prototype.set=function(t,e){var n=this.array_[t];return this.array_[t]=e,n},e.prototype.iterator=function(){return new Ct(this)},e.prototype.get=function(t){if(t<0||t>=this.size())throw new n;return this.array_[t]},e.prototype.isEmpty=function(){return 0===this.array_.length},e.prototype.size=function(){return this.array_.length},e.prototype.toArray=function(){for(var t=[],e=0,n=this.array_.length;e<n;e++)t.push(this.array_[e]);return t},e.prototype.remove=function(t){for(var e=!1,n=0,i=this.array_.length;n<i;n++)if(this.array_[n]===t){this.array_.splice(n,1),e=!0;break}return e},e}(xt),Ct=function(t){function e(e){t.call(this),this.arrayList_=e,this.position_=0}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.next=function(){if(this.position_===this.arrayList_.size())throw new i;return this.arrayList_.get(this.position_++)},e.prototype.hasNext=function(){return this.position_<this.arrayList_.size()},e.prototype.set=function(t){return this.arrayList_.set(this.position_-1,t)},e.prototype.remove=function(){this.arrayList_.remove(this.arrayList_.get(this.position_))},e}(Et),St=function(t){function e(){if(t.call(this),0===arguments.length);else if(1===arguments.length){var e=arguments[0];this.ensureCapacity(e.length),this.add(e,!0)}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this.ensureCapacity(n.length),this.add(n,i)}}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={coordArrayType:{configurable:!0}};return n.coordArrayType.get=function(){return new Array(0).fill(null)},e.prototype.getCoordinate=function(t){return this.get(t)},e.prototype.addAll=function(){if(2===arguments.length){for(var e=arguments[0],n=arguments[1],i=!1,r=e.iterator();r.hasNext();)this.add(r.next(),n),i=!0;return i}return t.prototype.addAll.apply(this,arguments)},e.prototype.clone=function(){for(var e=t.prototype.clone.call(this),n=0;n<this.size();n++)e.add(n,this.get(n).copy());return e},e.prototype.toCoordinateArray=function(){return this.toArray(e.coordArrayType)},e.prototype.add=function(){if(1===arguments.length){var e=arguments[0];t.prototype.add.call(this,e)}else if(2===arguments.length){if(arguments[0]instanceof Array&&"boolean"==typeof arguments[1]){var n=arguments[0],i=arguments[1];return this.add(n,i,!0),!0}if(arguments[0]instanceof C&&"boolean"==typeof arguments[1]){var r=arguments[0];if(!arguments[1]&&this.size()>=1){if(this.get(this.size()-1).equals2D(r))return null}t.prototype.add.call(this,r)}else if(arguments[0]instanceof Object&&"boolean"==typeof arguments[1]){var o=arguments[0],s=arguments[1];return this.add(o,s),!0}}else if(3===arguments.length){if("boolean"==typeof arguments[2]&&arguments[0]instanceof Array&&"boolean"==typeof arguments[1]){var a=arguments[0],u=arguments[1];if(arguments[2])for(var l=0;l<a.length;l++)this.add(a[l],u);else for(var c=a.length-1;c>=0;c--)this.add(a[c],u);return!0}if("boolean"==typeof arguments[2]&&Number.isInteger(arguments[0])&&arguments[1]instanceof C){var p=arguments[0],h=arguments[1];if(!arguments[2]){var f=this.size();if(f>0){if(p>0){if(this.get(p-1).equals2D(h))return null}if(p<f){if(this.get(p).equals2D(h))return null}}}t.prototype.add.call(this,p,h)}}else if(4===arguments.length){var g=arguments[0],d=arguments[1],y=arguments[2],_=arguments[3],m=1;y>_&&(m=-1);for(var v=y;v!==_;v+=m)this.add(g[v],d);return!0}},e.prototype.closeRing=function(){this.size()>0&&this.add(new C(this.get(0)),!1)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},Object.defineProperties(e,n),e}(Nt),Lt=function(){},bt={ForwardComparator:{configurable:!0},BidirectionalComparator:{configurable:!0},coordArrayType:{configurable:!0}};bt.ForwardComparator.get=function(){return wt},bt.BidirectionalComparator.get=function(){return Ot},bt.coordArrayType.get=function(){return new Array(0).fill(null)},Lt.prototype.interfaces_=function(){return[]},Lt.prototype.getClass=function(){return Lt},Lt.isRing=function(t){return!(t.length<4)&&!!t[0].equals2D(t[t.length-1])},Lt.ptNotInList=function(t,e){for(var n=0;n<t.length;n++){var i=t[n];if(Lt.indexOf(i,e)<0)return i}return null},Lt.scroll=function(t,e){var n=Lt.indexOf(e,t);if(n<0)return null;var i=new Array(t.length).fill(null);Y.arraycopy(t,n,i,0,t.length-n),Y.arraycopy(t,0,i,t.length-n,n),Y.arraycopy(i,0,t,0,t.length)},Lt.equals=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];if(t===e)return!0;if(null===t||null===e)return!1;if(t.length!==e.length)return!1;for(var n=0;n<t.length;n++)if(!t[n].equals(e[n]))return!1;return!0}if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2];if(i===r)return!0;if(null===i||null===r)return!1;if(i.length!==r.length)return!1;for(var s=0;s<i.length;s++)if(0!==o.compare(i[s],r[s]))return!1;return!0}},Lt.intersection=function(t,e){for(var n=new St,i=0;i<t.length;i++)e.intersects(t[i])&&n.add(t[i],!0);return n.toCoordinateArray()},Lt.hasRepeatedPoints=function(t){for(var e=1;e<t.length;e++)if(t[e-1].equals(t[e]))return!0;return!1},Lt.removeRepeatedPoints=function(t){if(!Lt.hasRepeatedPoints(t))return t;return new St(t,!1).toCoordinateArray()},Lt.reverse=function(t){for(var e=t.length-1,n=Math.trunc(e/2),i=0;i<=n;i++){var r=t[i];t[i]=t[e-i],t[e-i]=r}},Lt.removeNull=function(t){for(var e=0,n=0;n<t.length;n++)null!==t[n]&&e++;var i=new Array(e).fill(null);if(0===e)return i;for(var r=0,o=0;o<t.length;o++)null!==t[o]&&(i[r++]=t[o]);return i},Lt.copyDeep=function(){if(1===arguments.length){for(var t=arguments[0],e=new Array(t.length).fill(null),n=0;n<t.length;n++)e[n]=new C(t[n]);return e}if(5===arguments.length)for(var i=arguments[0],r=arguments[1],o=arguments[2],s=arguments[3],a=arguments[4],u=0;u<a;u++)o[s+u]=new C(i[r+u])},Lt.isEqualReversed=function(t,e){for(var n=0;n<t.length;n++){var i=t[n],r=e[t.length-n-1];if(0!==i.compareTo(r))return!1}return!0},Lt.envelope=function(t){for(var e=new j,n=0;n<t.length;n++)e.expandToInclude(t[n]);return e},Lt.toCoordinateArray=function(t){return t.toArray(Lt.coordArrayType)},Lt.atLeastNCoordinatesOrNothing=function(t,e){return e.length>=t?e:[]},Lt.indexOf=function(t,e){for(var n=0;n<e.length;n++)if(t.equals(e[n]))return n;return-1},Lt.increasingDirection=function(t){for(var e=0;e<Math.trunc(t.length/2);e++){var n=t.length-1-e,i=t[e].compareTo(t[n]);if(0!==i)return i}return 1},Lt.compare=function(t,e){for(var n=0;n<t.length&&n<e.length;){var i=t[n].compareTo(e[n]);if(0!==i)return i;n++}return n<e.length?-1:n<t.length?1:0},Lt.minCoordinate=function(t){for(var e=null,n=0;n<t.length;n++)(null===e||e.compareTo(t[n])>0)&&(e=t[n]);return e},Lt.extract=function(t,e,n){e=R.clamp(e,0,t.length);var i=(n=R.clamp(n,-1,t.length))-e+1;n<0&&(i=0),e>=t.length&&(i=0),n<e&&(i=0);var r=new Array(i).fill(null);if(0===i)return r;for(var o=0,s=e;s<=n;s++)r[o++]=t[s];return r},Object.defineProperties(Lt,bt);var wt=function(){};wt.prototype.compare=function(t,e){return Lt.compare(t,e)},wt.prototype.interfaces_=function(){return[N]},wt.prototype.getClass=function(){return wt};var Ot=function(){};Ot.prototype.compare=function(t,e){var n=t,i=e;if(n.length<i.length)return-1;if(n.length>i.length)return 1;if(0===n.length)return 0;var r=Lt.compare(n,i);return Lt.isEqualReversed(n,i)?0:r},Ot.prototype.OLDcompare=function(t,e){var n=t,i=e;if(n.length<i.length)return-1;if(n.length>i.length)return 1;if(0===n.length)return 0;for(var r=Lt.increasingDirection(n),o=Lt.increasingDirection(i),s=r>0?0:n.length-1,a=o>0?0:n.length-1,u=0;u<n.length;u++){var l=n[s].compareTo(i[a]);if(0!==l)return l;s+=r,a+=o}return 0},Ot.prototype.interfaces_=function(){return[N]},Ot.prototype.getClass=function(){return Ot};var Tt=function(){};Tt.prototype.get=function(){},Tt.prototype.put=function(){},Tt.prototype.size=function(){},Tt.prototype.values=function(){},Tt.prototype.entrySet=function(){};var Rt=function(t){function e(){t.apply(this,arguments)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e}(Tt);(r.prototype=new Error).name="OperationNotSupported",(o.prototype=new It).contains=function(){};var Pt=function(t){function e(){t.call(this),this.array_=[],arguments[0]instanceof It&&this.addAll(arguments[0])}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.contains=function(t){for(var e=0,n=this.array_.length;e<n;e++){if(this.array_[e]===t)return!0}return!1},e.prototype.add=function(t){return!this.contains(t)&&(this.array_.push(t),!0)},e.prototype.addAll=function(t){for(var e=t.iterator();e.hasNext();)this.add(e.next());return!0},e.prototype.remove=function(t){throw new Error},e.prototype.size=function(){return this.array_.length},e.prototype.isEmpty=function(){return 0===this.array_.length},e.prototype.toArray=function(){for(var t=[],e=0,n=this.array_.length;e<n;e++)t.push(this.array_[e]);return t},e.prototype.iterator=function(){return new Dt(this)},e}(o),Dt=function(t){function e(e){t.call(this),this.hashSet_=e,this.position_=0}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.next=function(){if(this.position_===this.hashSet_.size())throw new i;return this.hashSet_.array_[this.position_++]},e.prototype.hasNext=function(){return this.position_<this.hashSet_.size()},e.prototype.remove=function(){throw new r},e}(Et),Mt=0;(p.prototype=new Rt).get=function(t){for(var e=this.root_;null!==e;){var n=t.compareTo(e.key);if(n<0)e=e.left;else{if(!(n>0))return e.value;e=e.right}}return null},p.prototype.put=function(t,e){if(null===this.root_)return this.root_={key:t,value:e,left:null,right:null,parent:null,color:Mt,getValue:function(){return this.value},getKey:function(){return this.key}},this.size_=1,null;var n,i,r=this.root_;do{if(n=r,(i=t.compareTo(r.key))<0)r=r.left;else{if(!(i>0)){var o=r.value;return r.value=e,o}r=r.right}}while(null!==r);var s={key:t,left:null,right:null,value:e,parent:n,color:Mt,getValue:function(){return this.value},getKey:function(){return this.key}};return i<0?n.left=s:n.right=s,this.fixAfterInsertion(s),this.size_++,null},p.prototype.fixAfterInsertion=function(t){for(t.color=1;null!=t&&t!==this.root_&&1===t.parent.color;)if(a(t)===l(a(a(t)))){var e=c(a(a(t)));1===s(e)?(u(a(t),Mt),u(e,Mt),u(a(a(t)),1),t=a(a(t))):(t===c(a(t))&&(t=a(t),this.rotateLeft(t)),u(a(t),Mt),u(a(a(t)),1),this.rotateRight(a(a(t))))}else{var n=l(a(a(t)));1===s(n)?(u(a(t),Mt),u(n,Mt),u(a(a(t)),1),t=a(a(t))):(t===l(a(t))&&(t=a(t),this.rotateRight(t)),u(a(t),Mt),u(a(a(t)),1),this.rotateLeft(a(a(t))))}this.root_.color=Mt},p.prototype.values=function(){var t=new Nt,e=this.getFirstEntry();if(null!==e)for(t.add(e.value);null!==(e=p.successor(e));)t.add(e.value);return t},p.prototype.entrySet=function(){var t=new Pt,e=this.getFirstEntry();if(null!==e)for(t.add(e);null!==(e=p.successor(e));)t.add(e);return t},p.prototype.rotateLeft=function(t){if(null!=t){var e=t.right;t.right=e.left,null!=e.left&&(e.left.parent=t),e.parent=t.parent,null===t.parent?this.root_=e:t.parent.left===t?t.parent.left=e:t.parent.right=e,e.left=t,t.parent=e}},p.prototype.rotateRight=function(t){if(null!=t){var e=t.left;t.left=e.right,null!=e.right&&(e.right.parent=t),e.parent=t.parent,null===t.parent?this.root_=e:t.parent.right===t?t.parent.right=e:t.parent.left=e,e.right=t,t.parent=e}},p.prototype.getFirstEntry=function(){var t=this.root_;if(null!=t)for(;null!=t.left;)t=t.left;return t},p.successor=function(t){if(null===t)return null;if(null!==t.right){for(var e=t.right;null!==e.left;)e=e.left;return e}for(var n=t.parent,i=t;null!==n&&i===n.right;)i=n,n=n.parent;return n},p.prototype.size=function(){return this.size_};var At=function(){};At.prototype.interfaces_=function(){return[]},At.prototype.getClass=function(){return At},h.prototype=new o,(f.prototype=new h).contains=function(t){for(var e=0,n=this.array_.length;e<n;e++){if(0===this.array_[e].compareTo(t))return!0}return!1},f.prototype.add=function(t){if(this.contains(t))return!1;for(var e=0,n=this.array_.length;e<n;e++){if(1===this.array_[e].compareTo(t))return this.array_.splice(e,0,t),!0}return this.array_.push(t),!0},f.prototype.addAll=function(t){for(var e=t.iterator();e.hasNext();)this.add(e.next());return!0},f.prototype.remove=function(t){throw new r},f.prototype.size=function(){return this.array_.length},f.prototype.isEmpty=function(){return 0===this.array_.length},f.prototype.toArray=function(){for(var t=[],e=0,n=this.array_.length;e<n;e++)t.push(this.array_[e]);return t},f.prototype.iterator=function(){return new Ft(this)};var Ft=function(t){this.treeSet_=t,this.position_=0};Ft.prototype.next=function(){if(this.position_===this.treeSet_.size())throw new i;return this.treeSet_.array_[this.position_++]},Ft.prototype.hasNext=function(){return this.position_<this.treeSet_.size()},Ft.prototype.remove=function(){throw new r};var Gt=function(){};Gt.sort=function(){var t,e,n,i,r=arguments[0];if(1===arguments.length)i=function(t,e){return t.compareTo(e)},r.sort(i);else if(2===arguments.length)n=arguments[1],i=function(t,e){return n.compare(t,e)},r.sort(i);else if(3===arguments.length){(e=r.slice(arguments[1],arguments[2])).sort();var o=r.slice(0,arguments[1]).concat(e,r.slice(arguments[2],r.length));for(r.splice(0,r.length),t=0;t<o.length;t++)r.push(o[t])}else if(4===arguments.length)for(e=r.slice(arguments[1],arguments[2]),n=arguments[3],i=function(t,e){return n.compare(t,e)},e.sort(i),o=r.slice(0,arguments[1]).concat(e,r.slice(arguments[2],r.length)),r.splice(0,r.length),t=0;t<o.length;t++)r.push(o[t])},Gt.asList=function(t){for(var e=new Nt,n=0,i=t.length;n<i;n++)e.add(t[n]);return e};var qt=function(){},Bt={P:{configurable:!0},L:{configurable:!0},A:{configurable:!0},FALSE:{configurable:!0},TRUE:{configurable:!0},DONTCARE:{configurable:!0},SYM_FALSE:{configurable:!0},SYM_TRUE:{configurable:!0},SYM_DONTCARE:{configurable:!0},SYM_P:{configurable:!0},SYM_L:{configurable:!0},SYM_A:{configurable:!0}};Bt.P.get=function(){return 0},Bt.L.get=function(){return 1},Bt.A.get=function(){return 2},Bt.FALSE.get=function(){return-1},Bt.TRUE.get=function(){return-2},Bt.DONTCARE.get=function(){return-3},Bt.SYM_FALSE.get=function(){return"F"},Bt.SYM_TRUE.get=function(){return"T"},Bt.SYM_DONTCARE.get=function(){return"*"},Bt.SYM_P.get=function(){return"0"},Bt.SYM_L.get=function(){return"1"},Bt.SYM_A.get=function(){return"2"},qt.prototype.interfaces_=function(){return[]},qt.prototype.getClass=function(){return qt},qt.toDimensionSymbol=function(t){switch(t){case qt.FALSE:return qt.SYM_FALSE;case qt.TRUE:return qt.SYM_TRUE;case qt.DONTCARE:return qt.SYM_DONTCARE;case qt.P:return qt.SYM_P;case qt.L:return qt.SYM_L;case qt.A:return qt.SYM_A}throw new m("Unknown dimension value: "+t)},qt.toDimensionValue=function(t){switch(A.toUpperCase(t)){case qt.SYM_FALSE:return qt.FALSE;case qt.SYM_TRUE:return qt.TRUE;case qt.SYM_DONTCARE:return qt.DONTCARE;case qt.SYM_P:return qt.P;case qt.SYM_L:return qt.L;case qt.SYM_A:return qt.A}throw new m("Unknown dimension symbol: "+t)},Object.defineProperties(qt,Bt);var Vt=function(){};Vt.prototype.filter=function(t){},Vt.prototype.interfaces_=function(){return[]},Vt.prototype.getClass=function(){return Vt};var Ut=function(){};Ut.prototype.filter=function(t,e){},Ut.prototype.isDone=function(){},Ut.prototype.isGeometryChanged=function(){},Ut.prototype.interfaces_=function(){return[]},Ut.prototype.getClass=function(){return Ut};var zt=function(t){function e(e,n){if(t.call(this,n),this._geometries=e||[],t.hasNullElements(this._geometries))throw new m("geometries must not contain null elements")}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.computeEnvelopeInternal=function(){for(var t=new j,e=0;e<this._geometries.length;e++)t.expandToInclude(this._geometries[e].getEnvelopeInternal());return t},e.prototype.getGeometryN=function(t){return this._geometries[t]},e.prototype.getSortIndex=function(){return t.SORTINDEX_GEOMETRYCOLLECTION},e.prototype.getCoordinates=function(){for(var t=new Array(this.getNumPoints()).fill(null),e=-1,n=0;n<this._geometries.length;n++)for(var i=this._geometries[n].getCoordinates(),r=0;r<i.length;r++)t[++e]=i[r];return t},e.prototype.getArea=function(){for(var t=0,e=0;e<this._geometries.length;e++)t+=this._geometries[e].getArea();return t},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];if(!this.isEquivalentClass(e))return!1;var i=e;if(this._geometries.length!==i._geometries.length)return!1;for(var r=0;r<this._geometries.length;r++)if(!this._geometries[r].equalsExact(i._geometries[r],n))return!1;return!0}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.normalize=function(){for(var t=0;t<this._geometries.length;t++)this._geometries[t].normalize();Gt.sort(this._geometries)},e.prototype.getCoordinate=function(){return this.isEmpty()?null:this._geometries[0].getCoordinate()},e.prototype.getBoundaryDimension=function(){for(var t=qt.FALSE,e=0;e<this._geometries.length;e++)t=Math.max(t,this._geometries[e].getBoundaryDimension());return t},e.prototype.getDimension=function(){for(var t=qt.FALSE,e=0;e<this._geometries.length;e++)t=Math.max(t,this._geometries[e].getDimension());return t},e.prototype.getLength=function(){for(var t=0,e=0;e<this._geometries.length;e++)t+=this._geometries[e].getLength();return t},e.prototype.getNumPoints=function(){for(var t=0,e=0;e<this._geometries.length;e++)t+=this._geometries[e].getNumPoints();return t},e.prototype.getNumGeometries=function(){return this._geometries.length},e.prototype.reverse=function(){for(var t=this._geometries.length,e=new Array(t).fill(null),n=0;n<this._geometries.length;n++)e[n]=this._geometries[n].reverse();return this.getFactory().createGeometryCollection(e)},e.prototype.compareToSameClass=function(){if(1===arguments.length){var t=arguments[0],e=new f(Gt.asList(this._geometries)),n=new f(Gt.asList(t._geometries));return this.compare(e,n)}if(2===arguments.length){for(var i=arguments[0],r=arguments[1],o=i,s=this.getNumGeometries(),a=o.getNumGeometries(),u=0;u<s&&u<a;){var l=this.getGeometryN(u),c=o.getGeometryN(u),p=l.compareToSameClass(c,r);if(0!==p)return p;u++}return u<s?1:u<a?-1:0}},e.prototype.apply=function(){if(T(arguments[0],ft))for(var t=arguments[0],e=0;e<this._geometries.length;e++)this._geometries[e].apply(t);else if(T(arguments[0],Ut)){var n=arguments[0];if(0===this._geometries.length)return null;for(var i=0;i<this._geometries.length&&(this._geometries[i].apply(n),!n.isDone());i++);n.isGeometryChanged()&&this.geometryChanged()}else if(T(arguments[0],Vt)){var r=arguments[0];r.filter(this);for(var o=0;o<this._geometries.length;o++)this._geometries[o].apply(r)}else if(T(arguments[0],lt)){var s=arguments[0];s.filter(this);for(var a=0;a<this._geometries.length;a++)this._geometries[a].apply(s)}},e.prototype.getBoundary=function(){return this.checkNotGeometryCollection(this),et.shouldNeverReachHere(),null},e.prototype.clone=function(){var e=t.prototype.clone.call(this);e._geometries=new Array(this._geometries.length).fill(null);for(var n=0;n<this._geometries.length;n++)e._geometries[n]=this._geometries[n].clone();return e},e.prototype.getGeometryType=function(){return"GeometryCollection"},e.prototype.copy=function(){for(var t=new Array(this._geometries.length).fill(null),n=0;n<t.length;n++)t[n]=this._geometries[n].copy();return new e(t,this._factory)},e.prototype.isEmpty=function(){for(var t=0;t<this._geometries.length;t++)if(!this._geometries[t].isEmpty())return!1;return!0},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return-0x4f07bcb1f857d800},Object.defineProperties(e,n),e}(ct),Xt=function(t){function e(){t.apply(this,arguments)}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.getSortIndex=function(){return ct.SORTINDEX_MULTILINESTRING},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];return!!this.isEquivalentClass(e)&&t.prototype.equalsExact.call(this,e,n)}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.getBoundaryDimension=function(){return this.isClosed()?qt.FALSE:0},e.prototype.isClosed=function(){if(this.isEmpty())return!1;for(var t=0;t<this._geometries.length;t++)if(!this._geometries[t].isClosed())return!1;return!0},e.prototype.getDimension=function(){return 1},e.prototype.reverse=function(){for(var t=this._geometries.length,e=new Array(t).fill(null),n=0;n<this._geometries.length;n++)e[t-1-n]=this._geometries[n].reverse();return this.getFactory().createMultiLineString(e)},e.prototype.getBoundary=function(){return new Yt(this).getBoundary()},e.prototype.getGeometryType=function(){return"MultiLineString"},e.prototype.copy=function(){for(var t=new Array(this._geometries.length).fill(null),n=0;n<t.length;n++)t[n]=this._geometries[n].copy();return new e(t,this._factory)},e.prototype.interfaces_=function(){return[At]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return 0x7155d2ab4afa8000},Object.defineProperties(e,n),e}(zt),Yt=function(){if(this._geom=null,this._geomFact=null,this._bnRule=null,this._endpointMap=null,1===arguments.length){var t=arguments[0],e=gt.MOD2_BOUNDARY_RULE;this._geom=t,this._geomFact=t.getFactory(),this._bnRule=e}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this._geom=n,this._geomFact=n.getFactory(),this._bnRule=i}};Yt.prototype.boundaryMultiLineString=function(t){if(this._geom.isEmpty())return this.getEmptyMultiPoint();var e=this.computeBoundaryCoordinates(t);return 1===e.length?this._geomFact.createPoint(e[0]):this._geomFact.createMultiPointFromCoords(e)},Yt.prototype.getBoundary=function(){return this._geom instanceof Kt?this.boundaryLineString(this._geom):this._geom instanceof Xt?this.boundaryMultiLineString(this._geom):this._geom.getBoundary()},Yt.prototype.boundaryLineString=function(t){if(this._geom.isEmpty())return this.getEmptyMultiPoint();if(t.isClosed()){return this._bnRule.isInBoundary(2)?t.getStartPoint():this._geomFact.createMultiPoint()}return this._geomFact.createMultiPoint([t.getStartPoint(),t.getEndPoint()])},Yt.prototype.getEmptyMultiPoint=function(){return this._geomFact.createMultiPoint()},Yt.prototype.computeBoundaryCoordinates=function(t){var e=new Nt;this._endpointMap=new p;for(var n=0;n<t.getNumGeometries();n++){var i=t.getGeometryN(n);0!==i.getNumPoints()&&(this.addEndpoint(i.getCoordinateN(0)),this.addEndpoint(i.getCoordinateN(i.getNumPoints()-1)))}for(var r=this._endpointMap.entrySet().iterator();r.hasNext();){var o=r.next(),s=o.getValue().count;this._bnRule.isInBoundary(s)&&e.add(o.getKey())}return Lt.toCoordinateArray(e)},Yt.prototype.addEndpoint=function(t){var e=this._endpointMap.get(t);null===e&&(e=new kt,this._endpointMap.put(t,e)),e.count++},Yt.prototype.interfaces_=function(){return[]},Yt.prototype.getClass=function(){return Yt},Yt.getBoundary=function(){if(1===arguments.length){var t=arguments[0];return new Yt(t).getBoundary()}if(2===arguments.length){var e=arguments[0],n=arguments[1];return new Yt(e,n).getBoundary()}};var kt=function(){this.count=null};kt.prototype.interfaces_=function(){return[]},kt.prototype.getClass=function(){return kt};var jt=function(){},Ht={NEWLINE:{configurable:!0},SIMPLE_ORDINATE_FORMAT:{configurable:!0}};jt.prototype.interfaces_=function(){return[]},jt.prototype.getClass=function(){return jt},jt.chars=function(t,e){for(var n=new Array(e).fill(null),i=0;i<e;i++)n[i]=t;return String(n)},jt.getStackTrace=function(){if(1===arguments.length){var t=arguments[0],e=new function(){},n=new function(){}(e);return t.printStackTrace(n),e.toString()}if(2===arguments.length){for(var i=arguments[0],r=arguments[1],o="",s=new function(){}(new function(){}(jt.getStackTrace(i))),a=0;a<r;a++)try{o+=s.readLine()+jt.NEWLINE}catch(t){if(!(t instanceof g))throw t;et.shouldNeverReachHere()}return o}},jt.split=function(t,e){for(var n=e.length,i=new Nt,r=""+t,o=r.indexOf(e);o>=0;){var s=r.substring(0,o);i.add(s),o=(r=r.substring(o+n)).indexOf(e)}r.length>0&&i.add(r);for(var a=new Array(i.size()).fill(null),u=0;u<a.length;u++)a[u]=i.get(u);return a},jt.toString=function(){if(1===arguments.length){var t=arguments[0];return jt.SIMPLE_ORDINATE_FORMAT.format(t)}},jt.spaces=function(t){return jt.chars(" ",t)},Ht.NEWLINE.get=function(){return Y.getProperty("line.separator")},Ht.SIMPLE_ORDINATE_FORMAT.get=function(){return new function(){}("0.#")},Object.defineProperties(jt,Ht);var Wt=function(){};Wt.prototype.interfaces_=function(){return[]},Wt.prototype.getClass=function(){return Wt},Wt.copyCoord=function(t,e,n,i){for(var r=Math.min(t.getDimension(),n.getDimension()),o=0;o<r;o++)n.setOrdinate(i,o,t.getOrdinate(e,o))},Wt.isRing=function(t){var e=t.size();return 0===e||!(e<=3)&&(t.getOrdinate(0,V.X)===t.getOrdinate(e-1,V.X)&&t.getOrdinate(0,V.Y)===t.getOrdinate(e-1,V.Y))},Wt.isEqual=function(t,e){var n=t.size();if(n!==e.size())return!1;for(var i=Math.min(t.getDimension(),e.getDimension()),r=0;r<n;r++)for(var o=0;o<i;o++){var s=t.getOrdinate(r,o),a=e.getOrdinate(r,o);if(t.getOrdinate(r,o)!==e.getOrdinate(r,o)&&(!v.isNaN(s)||!v.isNaN(a)))return!1}return!0},Wt.extend=function(t,e,n){var i=t.create(n,e.getDimension()),r=e.size();if(Wt.copy(e,0,i,0,r),r>0)for(var o=r;o<n;o++)Wt.copy(e,r-1,i,o,1);return i},Wt.reverse=function(t){for(var e=t.size()-1,n=Math.trunc(e/2),i=0;i<=n;i++)Wt.swap(t,i,e-i)},Wt.swap=function(t,e,n){if(e===n)return null;for(var i=0;i<t.getDimension();i++){var r=t.getOrdinate(e,i);t.setOrdinate(e,i,t.getOrdinate(n,i)),t.setOrdinate(n,i,r)}},Wt.copy=function(t,e,n,i,r){for(var o=0;o<r;o++)Wt.copyCoord(t,e+o,n,i+o)},Wt.toString=function(){if(1===arguments.length){var t=arguments[0],e=t.size();if(0===e)return"()";var n=t.getDimension(),i=new D;i.append("(");for(var r=0;r<e;r++){r>0&&i.append(" ");for(var o=0;o<n;o++)o>0&&i.append(","),i.append(jt.toString(t.getOrdinate(r,o)))}return i.append(")"),i.toString()}},Wt.ensureValidRing=function(t,e){var n=e.size();if(0===n)return e;if(n<=3)return Wt.createClosedRing(t,e,4);return e.getOrdinate(0,V.X)===e.getOrdinate(n-1,V.X)&&e.getOrdinate(0,V.Y)===e.getOrdinate(n-1,V.Y)?e:Wt.createClosedRing(t,e,n+1)},Wt.createClosedRing=function(t,e,n){var i=t.create(n,e.getDimension()),r=e.size();Wt.copy(e,0,i,0,r);for(var o=r;o<n;o++)Wt.copy(e,0,i,o,1);return i};var Kt=function(t){function e(e,n){t.call(this,n),this._points=null,this.init(e)}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.computeEnvelopeInternal=function(){return this.isEmpty()?new j:this._points.expandEnvelope(new j)},e.prototype.isRing=function(){return this.isClosed()&&this.isSimple()},e.prototype.getSortIndex=function(){return t.SORTINDEX_LINESTRING},e.prototype.getCoordinates=function(){return this._points.toCoordinateArray()},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];if(!this.isEquivalentClass(e))return!1;var i=e;if(this._points.size()!==i._points.size())return!1;for(var r=0;r<this._points.size();r++)if(!this.equal(this._points.getCoordinate(r),i._points.getCoordinate(r),n))return!1;return!0}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.normalize=function(){for(var t=0;t<Math.trunc(this._points.size()/2);t++){var e=this._points.size()-1-t;if(!this._points.getCoordinate(t).equals(this._points.getCoordinate(e)))return this._points.getCoordinate(t).compareTo(this._points.getCoordinate(e))>0&&Wt.reverse(this._points),null}},e.prototype.getCoordinate=function(){return this.isEmpty()?null:this._points.getCoordinate(0)},e.prototype.getBoundaryDimension=function(){return this.isClosed()?qt.FALSE:0},e.prototype.isClosed=function(){return!this.isEmpty()&&this.getCoordinateN(0).equals2D(this.getCoordinateN(this.getNumPoints()-1))},e.prototype.getEndPoint=function(){return this.isEmpty()?null:this.getPointN(this.getNumPoints()-1)},e.prototype.getDimension=function(){return 1},e.prototype.getLength=function(){return at.computeLength(this._points)},e.prototype.getNumPoints=function(){return this._points.size()},e.prototype.reverse=function(){var t=this._points.copy();Wt.reverse(t);return this.getFactory().createLineString(t)},e.prototype.compareToSameClass=function(){if(1===arguments.length){for(var t=arguments[0],e=0,n=0;e<this._points.size()&&n<t._points.size();){var i=this._points.getCoordinate(e).compareTo(t._points.getCoordinate(n));if(0!==i)return i;e++,n++}return e<this._points.size()?1:n<t._points.size()?-1:0}if(2===arguments.length){var r=arguments[0];return arguments[1].compare(this._points,r._points)}},e.prototype.apply=function(){if(T(arguments[0],ft))for(var t=arguments[0],e=0;e<this._points.size();e++)t.filter(this._points.getCoordinate(e));else if(T(arguments[0],Ut)){var n=arguments[0];if(0===this._points.size())return null;for(var i=0;i<this._points.size()&&(n.filter(this._points,i),!n.isDone());i++);n.isGeometryChanged()&&this.geometryChanged()}else if(T(arguments[0],Vt)){arguments[0].filter(this)}else if(T(arguments[0],lt)){arguments[0].filter(this)}},e.prototype.getBoundary=function(){return new Yt(this).getBoundary()},e.prototype.isEquivalentClass=function(t){return t instanceof e},e.prototype.clone=function(){var e=t.prototype.clone.call(this);return e._points=this._points.clone(),e},e.prototype.getCoordinateN=function(t){return this._points.getCoordinate(t)},e.prototype.getGeometryType=function(){return"LineString"},e.prototype.copy=function(){return new e(this._points.copy(),this._factory)},e.prototype.getCoordinateSequence=function(){return this._points},e.prototype.isEmpty=function(){return 0===this._points.size()},e.prototype.init=function(t){if(null===t&&(t=this.getFactory().getCoordinateSequenceFactory().create([])),1===t.size())throw new m("Invalid number of points in LineString (found "+t.size()+" - must be 0 or >= 2)");this._points=t},e.prototype.isCoordinate=function(t){for(var e=0;e<this._points.size();e++)if(this._points.getCoordinate(e).equals(t))return!0;return!1},e.prototype.getStartPoint=function(){return this.isEmpty()?null:this.getPointN(0)},e.prototype.getPointN=function(t){return this.getFactory().createPoint(this._points.getCoordinate(t))},e.prototype.interfaces_=function(){return[At]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return 0x2b2b51ba435c8e00},Object.defineProperties(e,n),e}(ct),Jt=function(){};Jt.prototype.interfaces_=function(){return[]},Jt.prototype.getClass=function(){return Jt};var Qt=function(t){function e(e,n){t.call(this,n),this._coordinates=e||null,this.init(this._coordinates)}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.computeEnvelopeInternal=function(){if(this.isEmpty())return new j;var t=new j;return t.expandToInclude(this._coordinates.getX(0),this._coordinates.getY(0)),t},e.prototype.getSortIndex=function(){return t.SORTINDEX_POINT},e.prototype.getCoordinates=function(){return this.isEmpty()?[]:[this.getCoordinate()]},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];return!!this.isEquivalentClass(e)&&(!(!this.isEmpty()||!e.isEmpty())||this.isEmpty()===e.isEmpty()&&this.equal(e.getCoordinate(),this.getCoordinate(),n))}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.normalize=function(){},e.prototype.getCoordinate=function(){return 0!==this._coordinates.size()?this._coordinates.getCoordinate(0):null},e.prototype.getBoundaryDimension=function(){return qt.FALSE},e.prototype.getDimension=function(){return 0},e.prototype.getNumPoints=function(){return this.isEmpty()?0:1},e.prototype.reverse=function(){return this.copy()},e.prototype.getX=function(){if(null===this.getCoordinate())throw new Error("getX called on empty Point");return this.getCoordinate().x},e.prototype.compareToSameClass=function(){if(1===arguments.length){var t=arguments[0];return this.getCoordinate().compareTo(t.getCoordinate())}if(2===arguments.length){var e=arguments[0];return arguments[1].compare(this._coordinates,e._coordinates)}},e.prototype.apply=function(){if(T(arguments[0],ft)){var t=arguments[0];if(this.isEmpty())return null;t.filter(this.getCoordinate())}else if(T(arguments[0],Ut)){var e=arguments[0];if(this.isEmpty())return null;e.filter(this._coordinates,0),e.isGeometryChanged()&&this.geometryChanged()}else if(T(arguments[0],Vt)){arguments[0].filter(this)}else if(T(arguments[0],lt)){arguments[0].filter(this)}},e.prototype.getBoundary=function(){return this.getFactory().createGeometryCollection(null)},e.prototype.clone=function(){var e=t.prototype.clone.call(this);return e._coordinates=this._coordinates.clone(),e},e.prototype.getGeometryType=function(){return"Point"},e.prototype.copy=function(){return new e(this._coordinates.copy(),this._factory)},e.prototype.getCoordinateSequence=function(){return this._coordinates},e.prototype.getY=function(){if(null===this.getCoordinate())throw new Error("getY called on empty Point");return this.getCoordinate().y},e.prototype.isEmpty=function(){return 0===this._coordinates.size()},e.prototype.init=function(t){null===t&&(t=this.getFactory().getCoordinateSequenceFactory().create([])),et.isTrue(t.size()<=1),this._coordinates=t},e.prototype.isSimple=function(){return!0},e.prototype.interfaces_=function(){return[Jt]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return 0x44077bad161cbc00},Object.defineProperties(e,n),e}(ct),Zt=function(){};Zt.prototype.interfaces_=function(){return[]},Zt.prototype.getClass=function(){return Zt};var $t=function(t){function e(e,n,i){if(t.call(this,i),this._shell=null,this._holes=null,null===e&&(e=this.getFactory().createLinearRing()),null===n&&(n=[]),t.hasNullElements(n))throw new m("holes must not contain null elements");if(e.isEmpty()&&t.hasNonEmptyElements(n))throw new m("shell is empty but holes are not");this._shell=e,this._holes=n}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.computeEnvelopeInternal=function(){return this._shell.getEnvelopeInternal()},e.prototype.getSortIndex=function(){return t.SORTINDEX_POLYGON},e.prototype.getCoordinates=function(){if(this.isEmpty())return[];for(var t=new Array(this.getNumPoints()).fill(null),e=-1,n=this._shell.getCoordinates(),i=0;i<n.length;i++)t[++e]=n[i];for(var r=0;r<this._holes.length;r++)for(var o=this._holes[r].getCoordinates(),s=0;s<o.length;s++)t[++e]=o[s];return t},e.prototype.getArea=function(){var t=0;t+=Math.abs(at.signedArea(this._shell.getCoordinateSequence()));for(var e=0;e<this._holes.length;e++)t-=Math.abs(at.signedArea(this._holes[e].getCoordinateSequence()));return t},e.prototype.isRectangle=function(){if(0!==this.getNumInteriorRing())return!1;if(null===this._shell)return!1;if(5!==this._shell.getNumPoints())return!1;for(var t=this._shell.getCoordinateSequence(),e=this.getEnvelopeInternal(),n=0;n<5;n++){var i=t.getX(n);if(i!==e.getMinX()&&i!==e.getMaxX())return!1;var r=t.getY(n);if(r!==e.getMinY()&&r!==e.getMaxY())return!1}for(var o=t.getX(0),s=t.getY(0),a=1;a<=4;a++){var u=t.getX(a),l=t.getY(a);if(u!==o===(l!==s))return!1;o=u,s=l}return!0},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];if(!this.isEquivalentClass(e))return!1;var i=e,r=this._shell,o=i._shell;if(!r.equalsExact(o,n))return!1;if(this._holes.length!==i._holes.length)return!1;for(var s=0;s<this._holes.length;s++)if(!this._holes[s].equalsExact(i._holes[s],n))return!1;return!0}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.normalize=function(){if(0===arguments.length){this.normalize(this._shell,!0);for(var t=0;t<this._holes.length;t++)this.normalize(this._holes[t],!1);Gt.sort(this._holes)}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(e.isEmpty())return null;var i=new Array(e.getCoordinates().length-1).fill(null);Y.arraycopy(e.getCoordinates(),0,i,0,i.length);var r=Lt.minCoordinate(e.getCoordinates());Lt.scroll(i,r),Y.arraycopy(i,0,e.getCoordinates(),0,i.length),e.getCoordinates()[i.length]=i[0],at.isCCW(e.getCoordinates())===n&&Lt.reverse(e.getCoordinates())}},e.prototype.getCoordinate=function(){return this._shell.getCoordinate()},e.prototype.getNumInteriorRing=function(){return this._holes.length},e.prototype.getBoundaryDimension=function(){return 1},e.prototype.getDimension=function(){return 2},e.prototype.getLength=function(){var t=0;t+=this._shell.getLength();for(var e=0;e<this._holes.length;e++)t+=this._holes[e].getLength();return t},e.prototype.getNumPoints=function(){for(var t=this._shell.getNumPoints(),e=0;e<this._holes.length;e++)t+=this._holes[e].getNumPoints();return t},e.prototype.reverse=function(){var t=this.copy();t._shell=this._shell.copy().reverse(),t._holes=new Array(this._holes.length).fill(null);for(var e=0;e<this._holes.length;e++)t._holes[e]=this._holes[e].copy().reverse();return t},e.prototype.convexHull=function(){return this.getExteriorRing().convexHull()},e.prototype.compareToSameClass=function(){if(1===arguments.length){var t=arguments[0],e=this._shell,n=t._shell;return e.compareToSameClass(n)}if(2===arguments.length){var i=arguments[0],r=arguments[1],o=i,s=this._shell,a=o._shell,u=s.compareToSameClass(a,r);if(0!==u)return u;for(var l=this.getNumInteriorRing(),c=o.getNumInteriorRing(),p=0;p<l&&p<c;){var h=this.getInteriorRingN(p),f=o.getInteriorRingN(p),g=h.compareToSameClass(f,r);if(0!==g)return g;p++}return p<l?1:p<c?-1:0}},e.prototype.apply=function(t){if(T(t,ft)){this._shell.apply(t);for(var e=0;e<this._holes.length;e++)this._holes[e].apply(t)}else if(T(t,Ut)){if(this._shell.apply(t),!t.isDone())for(var n=0;n<this._holes.length&&(this._holes[n].apply(t),!t.isDone());n++);t.isGeometryChanged()&&this.geometryChanged()}else if(T(t,Vt))t.filter(this);else if(T(t,lt)){t.filter(this),this._shell.apply(t);for(var i=0;i<this._holes.length;i++)this._holes[i].apply(t)}},e.prototype.getBoundary=function(){if(this.isEmpty())return this.getFactory().createMultiLineString();var t=new Array(this._holes.length+1).fill(null);t[0]=this._shell;for(var e=0;e<this._holes.length;e++)t[e+1]=this._holes[e];return t.length<=1?this.getFactory().createLinearRing(t[0].getCoordinateSequence()):this.getFactory().createMultiLineString(t)},e.prototype.clone=function(){var e=t.prototype.clone.call(this);e._shell=this._shell.clone(),e._holes=new Array(this._holes.length).fill(null);for(var n=0;n<this._holes.length;n++)e._holes[n]=this._holes[n].clone();return e},e.prototype.getGeometryType=function(){return"Polygon"},e.prototype.copy=function(){for(var t=this._shell.copy(),n=new Array(this._holes.length).fill(null),i=0;i<n.length;i++)n[i]=this._holes[i].copy();return new e(t,n,this._factory)},e.prototype.getExteriorRing=function(){return this._shell},e.prototype.isEmpty=function(){return this._shell.isEmpty()},e.prototype.getInteriorRingN=function(t){return this._holes[t]},e.prototype.interfaces_=function(){return[Zt]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return-0x307ffefd8dc97200},Object.defineProperties(e,n),e}(ct),te=function(t){function e(){t.apply(this,arguments)}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.getSortIndex=function(){return ct.SORTINDEX_MULTIPOINT},e.prototype.isValid=function(){return!0},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];return!!this.isEquivalentClass(e)&&t.prototype.equalsExact.call(this,e,n)}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.getCoordinate=function(){if(1===arguments.length){var e=arguments[0];return this._geometries[e].getCoordinate()}return t.prototype.getCoordinate.apply(this,arguments)},e.prototype.getBoundaryDimension=function(){return qt.FALSE},e.prototype.getDimension=function(){return 0},e.prototype.getBoundary=function(){return this.getFactory().createGeometryCollection(null)},e.prototype.getGeometryType=function(){return"MultiPoint"},e.prototype.copy=function(){for(var t=new Array(this._geometries.length).fill(null),n=0;n<t.length;n++)t[n]=this._geometries[n].copy();return new e(t,this._factory)},e.prototype.interfaces_=function(){return[Jt]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return-0x6fb1ed4162e0fc00},Object.defineProperties(e,n),e}(zt),ee=function(t){function e(e,n){e instanceof C&&n instanceof _e&&(e=n.getCoordinateSequenceFactory().create(e)),t.call(this,e,n),this.validateConstruction()}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={MINIMUM_VALID_SIZE:{configurable:!0},serialVersionUID:{configurable:!0}};return e.prototype.getSortIndex=function(){return ct.SORTINDEX_LINEARRING},e.prototype.getBoundaryDimension=function(){return qt.FALSE},e.prototype.isClosed=function(){return!!this.isEmpty()||t.prototype.isClosed.call(this)},e.prototype.reverse=function(){var t=this._points.copy();Wt.reverse(t);return this.getFactory().createLinearRing(t)},e.prototype.validateConstruction=function(){if(!this.isEmpty()&&!t.prototype.isClosed.call(this))throw new m("Points of LinearRing do not form a closed linestring");if(this.getCoordinateSequence().size()>=1&&this.getCoordinateSequence().size()<e.MINIMUM_VALID_SIZE)throw new m("Invalid number of points in LinearRing (found "+this.getCoordinateSequence().size()+" - must be 0 or >= 4)")},e.prototype.getGeometryType=function(){return"LinearRing"},e.prototype.copy=function(){return new e(this._points.copy(),this._factory)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},n.MINIMUM_VALID_SIZE.get=function(){return 4},n.serialVersionUID.get=function(){return-0x3b229e262367a600},Object.defineProperties(e,n),e}(Kt),ne=function(t){function e(){t.apply(this,arguments)}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={serialVersionUID:{configurable:!0}};return e.prototype.getSortIndex=function(){return ct.SORTINDEX_MULTIPOLYGON},e.prototype.equalsExact=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];return!!this.isEquivalentClass(e)&&t.prototype.equalsExact.call(this,e,n)}return t.prototype.equalsExact.apply(this,arguments)},e.prototype.getBoundaryDimension=function(){return 1},e.prototype.getDimension=function(){return 2},e.prototype.reverse=function(){for(var t=this._geometries.length,e=new Array(t).fill(null),n=0;n<this._geometries.length;n++)e[n]=this._geometries[n].reverse();return this.getFactory().createMultiPolygon(e)},e.prototype.getBoundary=function(){if(this.isEmpty())return this.getFactory().createMultiLineString();for(var t=new Nt,e=0;e<this._geometries.length;e++)for(var n=this._geometries[e].getBoundary(),i=0;i<n.getNumGeometries();i++)t.add(n.getGeometryN(i));var r=new Array(t.size()).fill(null);return this.getFactory().createMultiLineString(t.toArray(r))},e.prototype.getGeometryType=function(){return"MultiPolygon"},e.prototype.copy=function(){for(var t=new Array(this._geometries.length).fill(null),n=0;n<t.length;n++)t[n]=this._geometries[n].copy();return new e(t,this._factory)},e.prototype.interfaces_=function(){return[Zt]},e.prototype.getClass=function(){return e},n.serialVersionUID.get=function(){return-0x7a5aa1369171980},Object.defineProperties(e,n),e}(zt),ie=function(t){this._factory=t||null,this._isUserDataCopied=!1},re={NoOpGeometryOperation:{configurable:!0},CoordinateOperation:{configurable:!0},CoordinateSequenceOperation:{configurable:!0}};ie.prototype.setCopyUserData=function(t){this._isUserDataCopied=t},ie.prototype.edit=function(t,e){if(null===t)return null;var n=this.editInternal(t,e);return this._isUserDataCopied&&n.setUserData(t.getUserData()),n},ie.prototype.editInternal=function(t,e){return null===this._factory&&(this._factory=t.getFactory()),t instanceof zt?this.editGeometryCollection(t,e):t instanceof $t?this.editPolygon(t,e):t instanceof Qt?e.edit(t,this._factory):t instanceof Kt?e.edit(t,this._factory):(et.shouldNeverReachHere("Unsupported Geometry class: "+t.getClass().getName()),null)},ie.prototype.editGeometryCollection=function(t,e){for(var n=e.edit(t,this._factory),i=new Nt,r=0;r<n.getNumGeometries();r++){var o=this.edit(n.getGeometryN(r),e);null===o||o.isEmpty()||i.add(o)}return n.getClass()===te?this._factory.createMultiPoint(i.toArray([])):n.getClass()===Xt?this._factory.createMultiLineString(i.toArray([])):n.getClass()===ne?this._factory.createMultiPolygon(i.toArray([])):this._factory.createGeometryCollection(i.toArray([]))},ie.prototype.editPolygon=function(t,e){var n=e.edit(t,this._factory);if(null===n&&(n=this._factory.createPolygon(null)),n.isEmpty())return n;var i=this.edit(n.getExteriorRing(),e);if(null===i||i.isEmpty())return this._factory.createPolygon();for(var r=new Nt,o=0;o<n.getNumInteriorRing();o++){var s=this.edit(n.getInteriorRingN(o),e);null===s||s.isEmpty()||r.add(s)}return this._factory.createPolygon(i,r.toArray([]))},ie.prototype.interfaces_=function(){return[]},ie.prototype.getClass=function(){return ie},ie.GeometryEditorOperation=function(){},re.NoOpGeometryOperation.get=function(){return oe},re.CoordinateOperation.get=function(){return se},re.CoordinateSequenceOperation.get=function(){return ae},Object.defineProperties(ie,re);var oe=function(){};oe.prototype.edit=function(t,e){return t},oe.prototype.interfaces_=function(){return[ie.GeometryEditorOperation]},oe.prototype.getClass=function(){return oe};var se=function(){};se.prototype.edit=function(t,e){var n=this.editCoordinates(t.getCoordinates(),t);return null===n?t:t instanceof ee?e.createLinearRing(n):t instanceof Kt?e.createLineString(n):t instanceof Qt?n.length>0?e.createPoint(n[0]):e.createPoint():t},se.prototype.interfaces_=function(){return[ie.GeometryEditorOperation]},se.prototype.getClass=function(){return se};var ae=function(){};ae.prototype.edit=function(t,e){return t instanceof ee?e.createLinearRing(this.edit(t.getCoordinateSequence(),t)):t instanceof Kt?e.createLineString(this.edit(t.getCoordinateSequence(),t)):t instanceof Qt?e.createPoint(this.edit(t.getCoordinateSequence(),t)):t},ae.prototype.interfaces_=function(){return[ie.GeometryEditorOperation]},ae.prototype.getClass=function(){return ae};var ue=function(){if(this._dimension=3,this._coordinates=null,1===arguments.length){if(arguments[0]instanceof Array)this._coordinates=arguments[0],this._dimension=3;else if(Number.isInteger(arguments[0])){var t=arguments[0];this._coordinates=new Array(t).fill(null);for(var e=0;e<t;e++)this._coordinates[e]=new C}else if(T(arguments[0],V)){var n=arguments[0];if(null===n)return this._coordinates=new Array(0).fill(null),null;this._dimension=n.getDimension(),this._coordinates=new Array(n.size()).fill(null);for(var i=0;i<this._coordinates.length;i++)this._coordinates[i]=n.getCoordinateCopy(i)}}else if(2===arguments.length)if(arguments[0]instanceof Array&&Number.isInteger(arguments[1])){var r=arguments[0],o=arguments[1];this._coordinates=r,this._dimension=o,null===r&&(this._coordinates=new Array(0).fill(null))}else if(Number.isInteger(arguments[0])&&Number.isInteger(arguments[1])){var s=arguments[0],a=arguments[1];this._coordinates=new Array(s).fill(null),this._dimension=a;for(var u=0;u<s;u++)this._coordinates[u]=new C}},le={serialVersionUID:{configurable:!0}};ue.prototype.setOrdinate=function(t,e,n){switch(e){case V.X:this._coordinates[t].x=n;break;case V.Y:this._coordinates[t].y=n;break;case V.Z:this._coordinates[t].z=n;break;default:throw new m("invalid ordinateIndex")}},ue.prototype.size=function(){return this._coordinates.length},ue.prototype.getOrdinate=function(t,e){switch(e){case V.X:return this._coordinates[t].x;case V.Y:return this._coordinates[t].y;case V.Z:return this._coordinates[t].z}return v.NaN},ue.prototype.getCoordinate=function(){if(1===arguments.length){var t=arguments[0];return this._coordinates[t]}if(2===arguments.length){var e=arguments[0],n=arguments[1];n.x=this._coordinates[e].x,n.y=this._coordinates[e].y,n.z=this._coordinates[e].z}},ue.prototype.getCoordinateCopy=function(t){return new C(this._coordinates[t])},ue.prototype.getDimension=function(){return this._dimension},ue.prototype.getX=function(t){return this._coordinates[t].x},ue.prototype.clone=function(){for(var t=new Array(this.size()).fill(null),e=0;e<this._coordinates.length;e++)t[e]=this._coordinates[e].clone();return new ue(t,this._dimension)},ue.prototype.expandEnvelope=function(t){for(var e=0;e<this._coordinates.length;e++)t.expandToInclude(this._coordinates[e]);return t},ue.prototype.copy=function(){for(var t=new Array(this.size()).fill(null),e=0;e<this._coordinates.length;e++)t[e]=this._coordinates[e].copy();return new ue(t,this._dimension)},ue.prototype.toString=function(){if(this._coordinates.length>0){var t=new D(17*this._coordinates.length);t.append("("),t.append(this._coordinates[0]);for(var e=1;e<this._coordinates.length;e++)t.append(", "),t.append(this._coordinates[e]);return t.append(")"),t.toString()}return"()"},ue.prototype.getY=function(t){return this._coordinates[t].y},ue.prototype.toCoordinateArray=function(){return this._coordinates},ue.prototype.interfaces_=function(){return[V,e]},ue.prototype.getClass=function(){return ue},le.serialVersionUID.get=function(){return-0xcb44a778db18e00},Object.defineProperties(ue,le);var ce=function(){},pe={serialVersionUID:{configurable:!0},instanceObject:{configurable:!0}};ce.prototype.readResolve=function(){return ce.instance()},ce.prototype.create=function(){if(1===arguments.length){if(arguments[0]instanceof Array){var t=arguments[0];return new ue(t)}if(T(arguments[0],V)){var e=arguments[0];return new ue(e)}}else if(2===arguments.length){var n=arguments[0],i=arguments[1];return i>3&&(i=3),i<2?new ue(n):new ue(n,i)}},ce.prototype.interfaces_=function(){return[b,e]},ce.prototype.getClass=function(){return ce},ce.instance=function(){return ce.instanceObject},pe.serialVersionUID.get=function(){return-0x38e49fa6cf6f2e00},pe.instanceObject.get=function(){return new ce},Object.defineProperties(ce,pe);var he=function(t){function e(){t.call(this),this.map_=new Map}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.get=function(t){return this.map_.get(t)||null},e.prototype.put=function(t,e){return this.map_.set(t,e),e},e.prototype.values=function(){for(var t=new Nt,e=this.map_.values(),n=e.next();!n.done;)t.add(n.value),n=e.next();return t},e.prototype.entrySet=function(){var t=new Pt;return this.map_.entries().forEach(function(e){return t.add(e)}),t},e.prototype.size=function(){return this.map_.size()},e}(Tt),fe=function t(){if(this._modelType=null,this._scale=null,0===arguments.length)this._modelType=t.FLOATING;else if(1===arguments.length)if(arguments[0]instanceof de){var e=arguments[0];this._modelType=e,e===t.FIXED&&this.setScale(1)}else if("number"==typeof arguments[0]){var n=arguments[0];this._modelType=t.FIXED,this.setScale(n)}else if(arguments[0]instanceof t){var i=arguments[0];this._modelType=i._modelType,this._scale=i._scale}},ge={serialVersionUID:{configurable:!0},maximumPreciseValue:{configurable:!0}};fe.prototype.equals=function(t){if(!(t instanceof fe))return!1;var e=t;return this._modelType===e._modelType&&this._scale===e._scale},fe.prototype.compareTo=function(t){var e=t,n=this.getMaximumSignificantDigits(),i=e.getMaximumSignificantDigits();return new M(n).compareTo(new M(i))},fe.prototype.getScale=function(){return this._scale},fe.prototype.isFloating=function(){return this._modelType===fe.FLOATING||this._modelType===fe.FLOATING_SINGLE},fe.prototype.getType=function(){return this._modelType},fe.prototype.toString=function(){var t="UNKNOWN";return this._modelType===fe.FLOATING?t="Floating":this._modelType===fe.FLOATING_SINGLE?t="Floating-Single":this._modelType===fe.FIXED&&(t="Fixed (Scale="+this.getScale()+")"),t},fe.prototype.makePrecise=function(){if("number"==typeof arguments[0]){var t=arguments[0];if(v.isNaN(t))return t;if(this._modelType===fe.FLOATING_SINGLE){return t}return this._modelType===fe.FIXED?Math.round(t*this._scale)/this._scale:t}if(arguments[0]instanceof C){var e=arguments[0];if(this._modelType===fe.FLOATING)return null;e.x=this.makePrecise(e.x),e.y=this.makePrecise(e.y)}},fe.prototype.getMaximumSignificantDigits=function(){var t=16;return this._modelType===fe.FLOATING?t=16:this._modelType===fe.FLOATING_SINGLE?t=6:this._modelType===fe.FIXED&&(t=1+Math.trunc(Math.ceil(Math.log(this.getScale())/Math.log(10)))),t},fe.prototype.setScale=function(t){this._scale=Math.abs(t)},fe.prototype.interfaces_=function(){return[e,E]},fe.prototype.getClass=function(){return fe},fe.mostPrecise=function(t,e){return t.compareTo(e)>=0?t:e},ge.serialVersionUID.get=function(){return 0x6bee6404e9a25c00},ge.maximumPreciseValue.get=function(){return 9007199254740992},Object.defineProperties(fe,ge);var de=function t(e){this._name=e||null,t.nameToTypeMap.put(e,this)},ye={serialVersionUID:{configurable:!0},nameToTypeMap:{configurable:!0}};de.prototype.readResolve=function(){return de.nameToTypeMap.get(this._name)},de.prototype.toString=function(){return this._name},de.prototype.interfaces_=function(){return[e]},de.prototype.getClass=function(){return de},ye.serialVersionUID.get=function(){return-552860263173159e4},ye.nameToTypeMap.get=function(){return new he},Object.defineProperties(de,ye),fe.Type=de,fe.FIXED=new de("FIXED"),fe.FLOATING=new de("FLOATING"),fe.FLOATING_SINGLE=new de("FLOATING SINGLE");var _e=function t(){this._precisionModel=new fe,this._SRID=0,this._coordinateSequenceFactory=t.getDefaultCoordinateSequenceFactory(),0===arguments.length||(1===arguments.length?T(arguments[0],b)?this._coordinateSequenceFactory=arguments[0]:arguments[0]instanceof fe&&(this._precisionModel=arguments[0]):2===arguments.length?(this._precisionModel=arguments[0],this._SRID=arguments[1]):3===arguments.length&&(this._precisionModel=arguments[0],this._SRID=arguments[1],this._coordinateSequenceFactory=arguments[2]))},me={serialVersionUID:{configurable:!0}};_e.prototype.toGeometry=function(t){return t.isNull()?this.createPoint(null):t.getMinX()===t.getMaxX()&&t.getMinY()===t.getMaxY()?this.createPoint(new C(t.getMinX(),t.getMinY())):t.getMinX()===t.getMaxX()||t.getMinY()===t.getMaxY()?this.createLineString([new C(t.getMinX(),t.getMinY()),new C(t.getMaxX(),t.getMaxY())]):this.createPolygon(this.createLinearRing([new C(t.getMinX(),t.getMinY()),new C(t.getMinX(),t.getMaxY()),new C(t.getMaxX(),t.getMaxY()),new C(t.getMaxX(),t.getMinY()),new C(t.getMinX(),t.getMinY())]),null)},_e.prototype.createLineString=function(t){return t?t instanceof Array?new Kt(this.getCoordinateSequenceFactory().create(t),this):T(t,V)?new Kt(t,this):void 0:new Kt(this.getCoordinateSequenceFactory().create([]),this)},_e.prototype.createMultiLineString=function(){if(0===arguments.length)return new Xt(null,this);if(1===arguments.length){var t=arguments[0];return new Xt(t,this)}},_e.prototype.buildGeometry=function(t){for(var e=null,n=!1,i=!1,r=t.iterator();r.hasNext();){var o=r.next(),s=o.getClass();null===e&&(e=s),s!==e&&(n=!0),o.isGeometryCollectionOrDerived()&&(i=!0)}if(null===e)return this.createGeometryCollection();if(n||i)return this.createGeometryCollection(_e.toGeometryArray(t));var a=t.iterator().next();if(t.size()>1){if(a instanceof $t)return this.createMultiPolygon(_e.toPolygonArray(t));if(a instanceof Kt)return this.createMultiLineString(_e.toLineStringArray(t));if(a instanceof Qt)return this.createMultiPoint(_e.toPointArray(t));et.shouldNeverReachHere("Unhandled class: "+a.getClass().getName())}return a},_e.prototype.createMultiPointFromCoords=function(t){return this.createMultiPoint(null!==t?this.getCoordinateSequenceFactory().create(t):null)},_e.prototype.createPoint=function(){if(0===arguments.length)return this.createPoint(this.getCoordinateSequenceFactory().create([]));if(1===arguments.length){if(arguments[0]instanceof C){var t=arguments[0];return this.createPoint(null!==t?this.getCoordinateSequenceFactory().create([t]):null)}if(T(arguments[0],V)){var e=arguments[0];return new Qt(e,this)}}},_e.prototype.getCoordinateSequenceFactory=function(){return this._coordinateSequenceFactory},_e.prototype.createPolygon=function(){if(0===arguments.length)return new $t(null,null,this);if(1===arguments.length){if(T(arguments[0],V)){var t=arguments[0];return this.createPolygon(this.createLinearRing(t))}if(arguments[0]instanceof Array){var e=arguments[0];return this.createPolygon(this.createLinearRing(e))}if(arguments[0]instanceof ee){var n=arguments[0];return this.createPolygon(n,null)}}else if(2===arguments.length){var i=arguments[0],r=arguments[1];return new $t(i,r,this)}},_e.prototype.getSRID=function(){return this._SRID},_e.prototype.createGeometryCollection=function(){if(0===arguments.length)return new zt(null,this);if(1===arguments.length){var t=arguments[0];return new zt(t,this)}},_e.prototype.createGeometry=function(t){return new ie(this).edit(t,{edit:function(){if(2===arguments.length){var t=arguments[0];return this._coordinateSequenceFactory.create(t)}}})},_e.prototype.getPrecisionModel=function(){return this._precisionModel},_e.prototype.createLinearRing=function(){if(0===arguments.length)return this.createLinearRing(this.getCoordinateSequenceFactory().create([]));if(1===arguments.length){if(arguments[0]instanceof Array){var t=arguments[0];return this.createLinearRing(null!==t?this.getCoordinateSequenceFactory().create(t):null)}if(T(arguments[0],V)){var e=arguments[0];return new ee(e,this)}}},_e.prototype.createMultiPolygon=function(){if(0===arguments.length)return new ne(null,this);if(1===arguments.length){var t=arguments[0];return new ne(t,this)}},_e.prototype.createMultiPoint=function(){if(0===arguments.length)return new te(null,this);if(1===arguments.length){if(arguments[0]instanceof Array){var t=arguments[0];return new te(t,this)}if(arguments[0]instanceof Array){var e=arguments[0];return this.createMultiPoint(null!==e?this.getCoordinateSequenceFactory().create(e):null)}if(T(arguments[0],V)){var n=arguments[0];if(null===n)return this.createMultiPoint(new Array(0).fill(null));for(var i=new Array(n.size()).fill(null),r=0;r<n.size();r++){var o=this.getCoordinateSequenceFactory().create(1,n.getDimension());Wt.copy(n,r,o,0,1),i[r]=this.createPoint(o)}return this.createMultiPoint(i)}}},_e.prototype.interfaces_=function(){return[e]},_e.prototype.getClass=function(){return _e},_e.toMultiPolygonArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toGeometryArray=function(t){if(null===t)return null;var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.getDefaultCoordinateSequenceFactory=function(){return ce.instance()},_e.toMultiLineStringArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toLineStringArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toMultiPointArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toLinearRingArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toPointArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.toPolygonArray=function(t){var e=new Array(t.size()).fill(null);return t.toArray(e)},_e.createPointFromInternalCoord=function(t,e){return e.getPrecisionModel().makePrecise(t),e.getFactory().createPoint(t)},me.serialVersionUID.get=function(){return-0x5ea75f2051eeb400},Object.defineProperties(_e,me);var ve=["Point","MultiPoint","LineString","MultiLineString","Polygon","MultiPolygon"],Ie=function(t){this.geometryFactory=t||new _e};Ie.prototype.read=function(t){var e,n=(e="string"==typeof t?JSON.parse(t):t).type;if(!Ee[n])throw new Error("Unknown GeoJSON type: "+e.type);return-1!==ve.indexOf(n)?Ee[n].apply(this,[e.coordinates]):"GeometryCollection"===n?Ee[n].apply(this,[e.geometries]):Ee[n].apply(this,[e])},Ie.prototype.write=function(t){var e=t.getGeometryType();if(!xe[e])throw new Error("Geometry is not supported");return xe[e].apply(this,[t])};var Ee={Feature:function(t){var e={};for(var n in t)e[n]=t[n];if(t.geometry){var i=t.geometry.type;if(!Ee[i])throw new Error("Unknown GeoJSON type: "+t.type);e.geometry=this.read(t.geometry)}return t.bbox&&(e.bbox=Ee.bbox.apply(this,[t.bbox])),e},FeatureCollection:function(t){var e={};if(t.features){e.features=[];for(var n=0;n<t.features.length;++n)e.features.push(this.read(t.features[n]))}return t.bbox&&(e.bbox=this.parse.bbox.apply(this,[t.bbox])),e},coordinates:function(t){for(var e=[],n=0;n<t.length;++n){var i=t[n];e.push(new C(i[0],i[1]))}return e},bbox:function(t){return this.geometryFactory.createLinearRing([new C(t[0],t[1]),new C(t[2],t[1]),new C(t[2],t[3]),new C(t[0],t[3]),new C(t[0],t[1])])},Point:function(t){var e=new C(t[0],t[1]);return this.geometryFactory.createPoint(e)},MultiPoint:function(t){for(var e=[],n=0;n<t.length;++n)e.push(Ee.Point.apply(this,[t[n]]));return this.geometryFactory.createMultiPoint(e)},LineString:function(t){var e=Ee.coordinates.apply(this,[t]);return this.geometryFactory.createLineString(e)},MultiLineString:function(t){for(var e=[],n=0;n<t.length;++n)e.push(Ee.LineString.apply(this,[t[n]]));return this.geometryFactory.createMultiLineString(e)},Polygon:function(t){for(var e=Ee.coordinates.apply(this,[t[0]]),n=this.geometryFactory.createLinearRing(e),i=[],r=1;r<t.length;++r){var o=t[r],s=Ee.coordinates.apply(this,[o]),a=this.geometryFactory.createLinearRing(s);i.push(a)}return this.geometryFactory.createPolygon(n,i)},MultiPolygon:function(t){for(var e=[],n=0;n<t.length;++n){var i=t[n];e.push(Ee.Polygon.apply(this,[i]))}return this.geometryFactory.createMultiPolygon(e)},GeometryCollection:function(t){for(var e=[],n=0;n<t.length;++n){var i=t[n];e.push(this.read(i))}return this.geometryFactory.createGeometryCollection(e)}},xe={coordinate:function(t){return[t.x,t.y]},Point:function(t){return{type:"Point",coordinates:xe.coordinate.apply(this,[t.getCoordinate()])}},MultiPoint:function(t){for(var e=[],n=0;n<t._geometries.length;++n){var i=t._geometries[n],r=xe.Point.apply(this,[i]);e.push(r.coordinates)}return{type:"MultiPoint",coordinates:e}},LineString:function(t){for(var e=[],n=t.getCoordinates(),i=0;i<n.length;++i){var r=n[i];e.push(xe.coordinate.apply(this,[r]))}return{type:"LineString",coordinates:e}},MultiLineString:function(t){for(var e=[],n=0;n<t._geometries.length;++n){var i=t._geometries[n],r=xe.LineString.apply(this,[i]);e.push(r.coordinates)}return{type:"MultiLineString",coordinates:e}},Polygon:function(t){var e=[],n=xe.LineString.apply(this,[t._shell]);e.push(n.coordinates);for(var i=0;i<t._holes.length;++i){var r=t._holes[i],o=xe.LineString.apply(this,[r]);e.push(o.coordinates)}return{type:"Polygon",coordinates:e}},MultiPolygon:function(t){for(var e=[],n=0;n<t._geometries.length;++n){var i=t._geometries[n],r=xe.Polygon.apply(this,[i]);e.push(r.coordinates)}return{type:"MultiPolygon",coordinates:e}},GeometryCollection:function(t){for(var e=[],n=0;n<t._geometries.length;++n){var i=t._geometries[n],r=i.getGeometryType();e.push(xe[r].apply(this,[i]))}return{type:"GeometryCollection",geometries:e}}},Ne=function(t){this.geometryFactory=t||new _e,this.precisionModel=this.geometryFactory.getPrecisionModel(),this.parser=new Ie(this.geometryFactory)};Ne.prototype.read=function(t){var e=this.parser.read(t);return this.precisionModel.getType()===fe.FIXED&&this.reducePrecision(e),e},Ne.prototype.reducePrecision=function(t){var e,n;if(t.coordinate)this.precisionModel.makePrecise(t.coordinate);else if(t.points)for(e=0,n=t.points.length;e<n;e++)this.precisionModel.makePrecise(t.points[e]);else if(t.geometries)for(e=0,n=t.geometries.length;e<n;e++)this.reducePrecision(t.geometries[e])};var Ce=function(){this.parser=new Ie(this.geometryFactory)};Ce.prototype.write=function(t){return this.parser.write(t)};var Se=function(){},Le={ON:{configurable:!0},LEFT:{configurable:!0},RIGHT:{configurable:!0}};Se.prototype.interfaces_=function(){return[]},Se.prototype.getClass=function(){return Se},Se.opposite=function(t){return t===Se.LEFT?Se.RIGHT:t===Se.RIGHT?Se.LEFT:t},Le.ON.get=function(){return 0},Le.LEFT.get=function(){return 1},Le.RIGHT.get=function(){return 2},Object.defineProperties(Se,Le),(d.prototype=new Error).name="EmptyStackException",(y.prototype=new xt).add=function(t){return this.array_.push(t),!0},y.prototype.get=function(t){if(t<0||t>=this.size())throw new Error;return this.array_[t]},y.prototype.push=function(t){return this.array_.push(t),t},y.prototype.pop=function(t){if(0===this.array_.length)throw new d;return this.array_.pop()},y.prototype.peek=function(){if(0===this.array_.length)throw new d;return this.array_[this.array_.length-1]},y.prototype.empty=function(){return 0===this.array_.length},y.prototype.isEmpty=function(){return this.empty()},y.prototype.search=function(t){return this.array_.indexOf(t)},y.prototype.size=function(){return this.array_.length},y.prototype.toArray=function(){for(var t=[],e=0,n=this.array_.length;e<n;e++)t.push(this.array_[e]);return t};var be=function(){this._minIndex=-1,this._minCoord=null,this._minDe=null,this._orientedDe=null};be.prototype.getCoordinate=function(){return this._minCoord},be.prototype.getRightmostSide=function(t,e){var n=this.getRightmostSideOfSegment(t,e);return n<0&&(n=this.getRightmostSideOfSegment(t,e-1)),n<0&&(this._minCoord=null,this.checkForRightmostCoordinate(t)),n},be.prototype.findRightmostEdgeAtVertex=function(){var t=this._minDe.getEdge().getCoordinates();et.isTrue(this._minIndex>0&&this._minIndex<t.length,"rightmost point expected to be interior vertex of edge");var e=t[this._minIndex-1],n=t[this._minIndex+1],i=at.computeOrientation(this._minCoord,n,e),r=!1;e.y<this._minCoord.y&&n.y<this._minCoord.y&&i===at.COUNTERCLOCKWISE?r=!0:e.y>this._minCoord.y&&n.y>this._minCoord.y&&i===at.CLOCKWISE&&(r=!0),r&&(this._minIndex=this._minIndex-1)},be.prototype.getRightmostSideOfSegment=function(t,e){var n=t.getEdge().getCoordinates();if(e<0||e+1>=n.length)return-1;if(n[e].y===n[e+1].y)return-1;var i=Se.LEFT;return n[e].y<n[e+1].y&&(i=Se.RIGHT),i},be.prototype.getEdge=function(){return this._orientedDe},be.prototype.checkForRightmostCoordinate=function(t){for(var e=t.getEdge().getCoordinates(),n=0;n<e.length-1;n++)(null===this._minCoord||e[n].x>this._minCoord.x)&&(this._minDe=t,this._minIndex=n,this._minCoord=e[n])},be.prototype.findRightmostEdgeAtNode=function(){var t=this._minDe.getNode().getEdges();this._minDe=t.getRightmostEdge(),this._minDe.isForward()||(this._minDe=this._minDe.getSym(),this._minIndex=this._minDe.getEdge().getCoordinates().length-1)},be.prototype.findEdge=function(t){for(var e=t.iterator();e.hasNext();){var n=e.next();n.isForward()&&this.checkForRightmostCoordinate(n)}et.isTrue(0!==this._minIndex||this._minCoord.equals(this._minDe.getCoordinate()),"inconsistency in rightmost processing"),0===this._minIndex?this.findRightmostEdgeAtNode():this.findRightmostEdgeAtVertex(),this._orientedDe=this._minDe;this.getRightmostSide(this._minDe,this._minIndex)===Se.LEFT&&(this._orientedDe=this._minDe.getSym())},be.prototype.interfaces_=function(){return[]},be.prototype.getClass=function(){return be};var we=function(t){function e(n,i){t.call(this,e.msgWithCoord(n,i)),this.pt=i?new C(i):null,this.name="TopologyException"}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.getCoordinate=function(){return this.pt},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e.msgWithCoord=function(t,e){return e?t:t+" [ "+e+" ]"},e}($),Oe=function(){this.array_=[]};Oe.prototype.addLast=function(t){this.array_.push(t)},Oe.prototype.removeFirst=function(){return this.array_.shift()},Oe.prototype.isEmpty=function(){return 0===this.array_.length};var Te=function(){this._finder=null,this._dirEdgeList=new Nt,this._nodes=new Nt,this._rightMostCoord=null,this._env=null,this._finder=new be};Te.prototype.clearVisitedEdges=function(){for(var t=this._dirEdgeList.iterator();t.hasNext();){t.next().setVisited(!1)}},Te.prototype.getRightmostCoordinate=function(){return this._rightMostCoord},Te.prototype.computeNodeDepth=function(t){for(var e=null,n=t.getEdges().iterator();n.hasNext();){var i=n.next();if(i.isVisited()||i.getSym().isVisited()){e=i;break}}if(null===e)throw new we("unable to find edge to compute depths at "+t.getCoordinate());t.getEdges().computeDepths(e);for(var r=t.getEdges().iterator();r.hasNext();){var o=r.next();o.setVisited(!0),this.copySymDepths(o)}},Te.prototype.computeDepth=function(t){this.clearVisitedEdges();var e=this._finder.getEdge();e.setEdgeDepths(Se.RIGHT,t),this.copySymDepths(e),this.computeDepths(e)},Te.prototype.create=function(t){this.addReachable(t),this._finder.findEdge(this._dirEdgeList),this._rightMostCoord=this._finder.getCoordinate()},Te.prototype.findResultEdges=function(){for(var t=this._dirEdgeList.iterator();t.hasNext();){var e=t.next();e.getDepth(Se.RIGHT)>=1&&e.getDepth(Se.LEFT)<=0&&!e.isInteriorAreaEdge()&&e.setInResult(!0)}},Te.prototype.computeDepths=function(t){var e=new Pt,n=new Oe,i=t.getNode();for(n.addLast(i),e.add(i),t.setVisited(!0);!n.isEmpty();){var r=n.removeFirst();e.add(r),this.computeNodeDepth(r);for(var o=r.getEdges().iterator();o.hasNext();){var s=o.next().getSym();if(!s.isVisited()){var a=s.getNode();e.contains(a)||(n.addLast(a),e.add(a))}}}},Te.prototype.compareTo=function(t){var e=t;return this._rightMostCoord.x<e._rightMostCoord.x?-1:this._rightMostCoord.x>e._rightMostCoord.x?1:0},Te.prototype.getEnvelope=function(){if(null===this._env){for(var t=new j,e=this._dirEdgeList.iterator();e.hasNext();)for(var n=e.next().getEdge().getCoordinates(),i=0;i<n.length-1;i++)t.expandToInclude(n[i]);this._env=t}return this._env},Te.prototype.addReachable=function(t){var e=new y;for(e.add(t);!e.empty();){var n=e.pop();this.add(n,e)}},Te.prototype.copySymDepths=function(t){var e=t.getSym();e.setDepth(Se.LEFT,t.getDepth(Se.RIGHT)),e.setDepth(Se.RIGHT,t.getDepth(Se.LEFT))},Te.prototype.add=function(t,e){t.setVisited(!0),this._nodes.add(t);for(var n=t.getEdges().iterator();n.hasNext();){var i=n.next();this._dirEdgeList.add(i);var r=i.getSym().getNode();r.isVisited()||e.push(r)}},Te.prototype.getNodes=function(){return this._nodes},Te.prototype.getDirectedEdges=function(){return this._dirEdgeList},Te.prototype.interfaces_=function(){return[E]},Te.prototype.getClass=function(){return Te};var Re=function t(){if(this.location=null,1===arguments.length){if(arguments[0]instanceof Array){var e=arguments[0];this.init(e.length)}else if(Number.isInteger(arguments[0])){var n=arguments[0];this.init(1),this.location[Se.ON]=n}else if(arguments[0]instanceof t){var i=arguments[0];if(this.init(i.location.length),null!==i)for(var r=0;r<this.location.length;r++)this.location[r]=i.location[r]}}else if(3===arguments.length){var o=arguments[0],s=arguments[1],a=arguments[2];this.init(3),this.location[Se.ON]=o,this.location[Se.LEFT]=s,this.location[Se.RIGHT]=a}};Re.prototype.setAllLocations=function(t){for(var e=0;e<this.location.length;e++)this.location[e]=t},Re.prototype.isNull=function(){for(var t=0;t<this.location.length;t++)if(this.location[t]!==w.NONE)return!1;return!0},Re.prototype.setAllLocationsIfNull=function(t){for(var e=0;e<this.location.length;e++)this.location[e]===w.NONE&&(this.location[e]=t)},Re.prototype.isLine=function(){return 1===this.location.length},Re.prototype.merge=function(t){if(t.location.length>this.location.length){var e=new Array(3).fill(null);e[Se.ON]=this.location[Se.ON],e[Se.LEFT]=w.NONE,e[Se.RIGHT]=w.NONE,this.location=e}for(var n=0;n<this.location.length;n++)this.location[n]===w.NONE&&n<t.location.length&&(this.location[n]=t.location[n])},Re.prototype.getLocations=function(){return this.location},Re.prototype.flip=function(){if(this.location.length<=1)return null;var t=this.location[Se.LEFT];this.location[Se.LEFT]=this.location[Se.RIGHT],this.location[Se.RIGHT]=t},Re.prototype.toString=function(){var t=new D;return this.location.length>1&&t.append(w.toLocationSymbol(this.location[Se.LEFT])),t.append(w.toLocationSymbol(this.location[Se.ON])),this.location.length>1&&t.append(w.toLocationSymbol(this.location[Se.RIGHT])),t.toString()},Re.prototype.setLocations=function(t,e,n){this.location[Se.ON]=t,this.location[Se.LEFT]=e,this.location[Se.RIGHT]=n},Re.prototype.get=function(t){return t<this.location.length?this.location[t]:w.NONE},Re.prototype.isArea=function(){return this.location.length>1},Re.prototype.isAnyNull=function(){for(var t=0;t<this.location.length;t++)if(this.location[t]===w.NONE)return!0;return!1},Re.prototype.setLocation=function(){if(1===arguments.length){var t=arguments[0];this.setLocation(Se.ON,t)}else if(2===arguments.length){var e=arguments[0],n=arguments[1];this.location[e]=n}},Re.prototype.init=function(t){this.location=new Array(t).fill(null),this.setAllLocations(w.NONE)},Re.prototype.isEqualOnSide=function(t,e){return this.location[e]===t.location[e]},Re.prototype.allPositionsEqual=function(t){for(var e=0;e<this.location.length;e++)if(this.location[e]!==t)return!1;return!0},Re.prototype.interfaces_=function(){return[]},Re.prototype.getClass=function(){return Re};var Pe=function t(){if(this.elt=new Array(2).fill(null),1===arguments.length){if(Number.isInteger(arguments[0])){var e=arguments[0];this.elt[0]=new Re(e),this.elt[1]=new Re(e)}else if(arguments[0]instanceof t){var n=arguments[0];this.elt[0]=new Re(n.elt[0]),this.elt[1]=new Re(n.elt[1])}}else if(2===arguments.length){var i=arguments[0],r=arguments[1];this.elt[0]=new Re(w.NONE),this.elt[1]=new Re(w.NONE),this.elt[i].setLocation(r)}else if(3===arguments.length){var o=arguments[0],s=arguments[1],a=arguments[2];this.elt[0]=new Re(o,s,a),this.elt[1]=new Re(o,s,a)}else if(4===arguments.length){var u=arguments[0],l=arguments[1],c=arguments[2],p=arguments[3];this.elt[0]=new Re(w.NONE,w.NONE,w.NONE),this.elt[1]=new Re(w.NONE,w.NONE,w.NONE),this.elt[u].setLocations(l,c,p)}};Pe.prototype.getGeometryCount=function(){var t=0;return this.elt[0].isNull()||t++,this.elt[1].isNull()||t++,t},Pe.prototype.setAllLocations=function(t,e){this.elt[t].setAllLocations(e)},Pe.prototype.isNull=function(t){return this.elt[t].isNull()},Pe.prototype.setAllLocationsIfNull=function(){if(1===arguments.length){var t=arguments[0];this.setAllLocationsIfNull(0,t),this.setAllLocationsIfNull(1,t)}else if(2===arguments.length){var e=arguments[0],n=arguments[1];this.elt[e].setAllLocationsIfNull(n)}},Pe.prototype.isLine=function(t){return this.elt[t].isLine()},Pe.prototype.merge=function(t){for(var e=0;e<2;e++)null===this.elt[e]&&null!==t.elt[e]?this.elt[e]=new Re(t.elt[e]):this.elt[e].merge(t.elt[e])},Pe.prototype.flip=function(){this.elt[0].flip(),this.elt[1].flip()},Pe.prototype.getLocation=function(){if(1===arguments.length){var t=arguments[0];return this.elt[t].get(Se.ON)}if(2===arguments.length){var e=arguments[0],n=arguments[1];return this.elt[e].get(n)}},Pe.prototype.toString=function(){var t=new D;return null!==this.elt[0]&&(t.append("A:"),t.append(this.elt[0].toString())),null!==this.elt[1]&&(t.append(" B:"),t.append(this.elt[1].toString())),t.toString()},Pe.prototype.isArea=function(){if(0===arguments.length)return this.elt[0].isArea()||this.elt[1].isArea();if(1===arguments.length){var t=arguments[0];return this.elt[t].isArea()}},Pe.prototype.isAnyNull=function(t){return this.elt[t].isAnyNull()},Pe.prototype.setLocation=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];this.elt[t].setLocation(Se.ON,e)}else if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2];this.elt[n].setLocation(i,r)}},Pe.prototype.isEqualOnSide=function(t,e){return this.elt[0].isEqualOnSide(t.elt[0],e)&&this.elt[1].isEqualOnSide(t.elt[1],e)},Pe.prototype.allPositionsEqual=function(t,e){return this.elt[t].allPositionsEqual(e)},Pe.prototype.toLine=function(t){this.elt[t].isArea()&&(this.elt[t]=new Re(this.elt[t].location[0]))},Pe.prototype.interfaces_=function(){return[]},Pe.prototype.getClass=function(){return Pe},Pe.toLineLabel=function(t){for(var e=new Pe(w.NONE),n=0;n<2;n++)e.setLocation(n,t.getLocation(n));return e};var De=function(){this._startDe=null,this._maxNodeDegree=-1,this._edges=new Nt,this._pts=new Nt,this._label=new Pe(w.NONE),this._ring=null,this._isHole=null,this._shell=null,this._holes=new Nt,this._geometryFactory=null;var t=arguments[0],e=arguments[1];this._geometryFactory=e,this.computePoints(t),this.computeRing()};De.prototype.computeRing=function(){if(null!==this._ring)return null;for(var t=new Array(this._pts.size()).fill(null),e=0;e<this._pts.size();e++)t[e]=this._pts.get(e);this._ring=this._geometryFactory.createLinearRing(t),this._isHole=at.isCCW(this._ring.getCoordinates())},De.prototype.isIsolated=function(){return 1===this._label.getGeometryCount()},De.prototype.computePoints=function(t){this._startDe=t;var e=t,n=!0;do{if(null===e)throw new we("Found null DirectedEdge");if(e.getEdgeRing()===this)throw new we("Directed Edge visited twice during ring-building at "+e.getCoordinate());this._edges.add(e);var i=e.getLabel();et.isTrue(i.isArea()),this.mergeLabel(i),this.addPoints(e.getEdge(),e.isForward(),n),n=!1,this.setEdgeRing(e,this),e=this.getNext(e)}while(e!==this._startDe)},De.prototype.getLinearRing=function(){return this._ring},De.prototype.getCoordinate=function(t){return this._pts.get(t)},De.prototype.computeMaxNodeDegree=function(){this._maxNodeDegree=0;var t=this._startDe;do{var e=t.getNode().getEdges().getOutgoingDegree(this);e>this._maxNodeDegree&&(this._maxNodeDegree=e),t=this.getNext(t)}while(t!==this._startDe);this._maxNodeDegree*=2},De.prototype.addPoints=function(t,e,n){var i=t.getCoordinates();if(e){var r=1;n&&(r=0);for(var o=r;o<i.length;o++)this._pts.add(i[o])}else{var s=i.length-2;n&&(s=i.length-1);for(var a=s;a>=0;a--)this._pts.add(i[a])}},De.prototype.isHole=function(){return this._isHole},De.prototype.setInResult=function(){var t=this._startDe;do{t.getEdge().setInResult(!0),t=t.getNext()}while(t!==this._startDe)},De.prototype.containsPoint=function(t){var e=this.getLinearRing();if(!e.getEnvelopeInternal().contains(t))return!1;if(!at.isPointInRing(t,e.getCoordinates()))return!1;for(var n=this._holes.iterator();n.hasNext();){if(n.next().containsPoint(t))return!1}return!0},De.prototype.addHole=function(t){this._holes.add(t)},De.prototype.isShell=function(){return null===this._shell},De.prototype.getLabel=function(){return this._label},De.prototype.getEdges=function(){return this._edges},De.prototype.getMaxNodeDegree=function(){return this._maxNodeDegree<0&&this.computeMaxNodeDegree(),this._maxNodeDegree},De.prototype.getShell=function(){return this._shell},De.prototype.mergeLabel=function(){if(1===arguments.length){var t=arguments[0];this.mergeLabel(t,0),this.mergeLabel(t,1)}else if(2===arguments.length){var e=arguments[0],n=arguments[1],i=e.getLocation(n,Se.RIGHT);if(i===w.NONE)return null;if(this._label.getLocation(n)===w.NONE)return this._label.setLocation(n,i),null}},De.prototype.setShell=function(t){this._shell=t,null!==t&&t.addHole(this)},De.prototype.toPolygon=function(t){for(var e=new Array(this._holes.size()).fill(null),n=0;n<this._holes.size();n++)e[n]=this._holes.get(n).getLinearRing();return t.createPolygon(this.getLinearRing(),e)},De.prototype.interfaces_=function(){return[]},De.prototype.getClass=function(){return De};var Me=function(t){function e(){var e=arguments[0],n=arguments[1];t.call(this,e,n)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.setEdgeRing=function(t,e){t.setMinEdgeRing(e)},e.prototype.getNext=function(t){return t.getNextMin()},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(De),Ae=function(t){function e(){var e=arguments[0],n=arguments[1];t.call(this,e,n)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.buildMinimalRings=function(){var t=new Nt,e=this._startDe;do{if(null===e.getMinEdgeRing()){var n=new Me(e,this._geometryFactory);t.add(n)}e=e.getNext()}while(e!==this._startDe);return t},e.prototype.setEdgeRing=function(t,e){t.setEdgeRing(e)},e.prototype.linkDirectedEdgesForMinimalEdgeRings=function(){var t=this._startDe;do{t.getNode().getEdges().linkMinimalDirectedEdges(this),t=t.getNext()}while(t!==this._startDe)},e.prototype.getNext=function(t){return t.getNext()},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(De),Fe=function(){if(this._label=null,this._isInResult=!1,this._isCovered=!1,this._isCoveredSet=!1,this._isVisited=!1,0===arguments.length);else if(1===arguments.length){var t=arguments[0];this._label=t}};Fe.prototype.setVisited=function(t){this._isVisited=t},Fe.prototype.setInResult=function(t){this._isInResult=t},Fe.prototype.isCovered=function(){return this._isCovered},Fe.prototype.isCoveredSet=function(){return this._isCoveredSet},Fe.prototype.setLabel=function(t){this._label=t},Fe.prototype.getLabel=function(){return this._label},Fe.prototype.setCovered=function(t){this._isCovered=t,this._isCoveredSet=!0},Fe.prototype.updateIM=function(t){et.isTrue(this._label.getGeometryCount()>=2,"found partial label"),this.computeIM(t)},Fe.prototype.isInResult=function(){return this._isInResult},Fe.prototype.isVisited=function(){return this._isVisited},Fe.prototype.interfaces_=function(){return[]},Fe.prototype.getClass=function(){return Fe};var Ge=function(t){function e(){t.call(this),this._coord=null,this._edges=null;var e=arguments[0],n=arguments[1];this._coord=e,this._edges=n,this._label=new Pe(0,w.NONE)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.isIncidentEdgeInResult=function(){for(var t=this.getEdges().getEdges().iterator();t.hasNext();){if(t.next().getEdge().isInResult())return!0}return!1},e.prototype.isIsolated=function(){return 1===this._label.getGeometryCount()},e.prototype.getCoordinate=function(){return this._coord},e.prototype.print=function(t){t.println("node "+this._coord+" lbl: "+this._label)},e.prototype.computeIM=function(t){},e.prototype.computeMergedLocation=function(t,e){var n=w.NONE;if(n=this._label.getLocation(e),!t.isNull(e)){var i=t.getLocation(e);n!==w.BOUNDARY&&(n=i)}return n},e.prototype.setLabel=function(){if(2!==arguments.length)return t.prototype.setLabel.apply(this,arguments);var e=arguments[0],n=arguments[1];null===this._label?this._label=new Pe(e,n):this._label.setLocation(e,n)},e.prototype.getEdges=function(){return this._edges},e.prototype.mergeLabel=function(){if(arguments[0]instanceof e){var t=arguments[0];this.mergeLabel(t._label)}else if(arguments[0]instanceof Pe)for(var n=arguments[0],i=0;i<2;i++){var r=this.computeMergedLocation(n,i);this._label.getLocation(i)===w.NONE&&this._label.setLocation(i,r)}},e.prototype.add=function(t){this._edges.insert(t),t.setNode(this)},e.prototype.setLabelBoundary=function(t){if(null===this._label)return null;var e=w.NONE;null!==this._label&&(e=this._label.getLocation(t));var n=null;switch(e){case w.BOUNDARY:n=w.INTERIOR;break;case w.INTERIOR:default:n=w.BOUNDARY}this._label.setLocation(t,n)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Fe),qe=function(){this.nodeMap=new p,this.nodeFact=null;var t=arguments[0];this.nodeFact=t};qe.prototype.find=function(t){return this.nodeMap.get(t)},qe.prototype.addNode=function(){if(arguments[0]instanceof C){var t=arguments[0],e=this.nodeMap.get(t);return null===e&&(e=this.nodeFact.createNode(t),this.nodeMap.put(t,e)),e}if(arguments[0]instanceof Ge){var n=arguments[0],i=this.nodeMap.get(n.getCoordinate());return null===i?(this.nodeMap.put(n.getCoordinate(),n),n):(i.mergeLabel(n),i)}},qe.prototype.print=function(t){for(var e=this.iterator();e.hasNext();){e.next().print(t)}},qe.prototype.iterator=function(){return this.nodeMap.values().iterator()},qe.prototype.values=function(){return this.nodeMap.values()},qe.prototype.getBoundaryNodes=function(t){for(var e=new Nt,n=this.iterator();n.hasNext();){var i=n.next();i.getLabel().getLocation(t)===w.BOUNDARY&&e.add(i)}return e},qe.prototype.add=function(t){var e=t.getCoordinate();this.addNode(e).add(t)},qe.prototype.interfaces_=function(){return[]},qe.prototype.getClass=function(){return qe};var Be=function(){},Ve={NE:{configurable:!0},NW:{configurable:!0},SW:{configurable:!0},SE:{configurable:!0}};Be.prototype.interfaces_=function(){return[]},Be.prototype.getClass=function(){return Be},Be.isNorthern=function(t){return t===Be.NE||t===Be.NW},Be.isOpposite=function(t,e){if(t===e)return!1;return 2===(t-e+4)%4},Be.commonHalfPlane=function(t,e){if(t===e)return t;if(2===(t-e+4)%4)return-1;var n=t<e?t:e;return 0===n&&3===(t>e?t:e)?3:n},Be.isInHalfPlane=function(t,e){return e===Be.SE?t===Be.SE||t===Be.SW:t===e||t===e+1},Be.quadrant=function(){if("number"==typeof arguments[0]&&"number"==typeof arguments[1]){var t=arguments[0],e=arguments[1];if(0===t&&0===e)throw new m("Cannot compute the quadrant for point ( "+t+", "+e+" )");return t>=0?e>=0?Be.NE:Be.SE:e>=0?Be.NW:Be.SW}if(arguments[0]instanceof C&&arguments[1]instanceof C){var n=arguments[0],i=arguments[1];if(i.x===n.x&&i.y===n.y)throw new m("Cannot compute the quadrant for two identical points "+n);return i.x>=n.x?i.y>=n.y?Be.NE:Be.SE:i.y>=n.y?Be.NW:Be.SW}},Ve.NE.get=function(){return 0},Ve.NW.get=function(){return 1},Ve.SW.get=function(){return 2},Ve.SE.get=function(){return 3},Object.defineProperties(Be,Ve);var Ue=function(){if(this._edge=null,this._label=null,this._node=null,this._p0=null,this._p1=null,this._dx=null,this._dy=null,this._quadrant=null,1===arguments.length){var t=arguments[0];this._edge=t}else if(3===arguments.length){var e=arguments[0],n=arguments[1],i=arguments[2];this._edge=e,this.init(n,i),this._label=null}else if(4===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2],a=arguments[3];this._edge=r,this.init(o,s),this._label=a}};Ue.prototype.compareDirection=function(t){return this._dx===t._dx&&this._dy===t._dy?0:this._quadrant>t._quadrant?1:this._quadrant<t._quadrant?-1:at.computeOrientation(t._p0,t._p1,this._p1)},Ue.prototype.getDy=function(){return this._dy},Ue.prototype.getCoordinate=function(){return this._p0},Ue.prototype.setNode=function(t){this._node=t},Ue.prototype.print=function(t){var e=Math.atan2(this._dy,this._dx),n=this.getClass().getName(),i=n.lastIndexOf("."),r=n.substring(i+1);t.print("  "+r+": "+this._p0+" - "+this._p1+" "+this._quadrant+":"+e+"   "+this._label)},Ue.prototype.compareTo=function(t){var e=t;return this.compareDirection(e)},Ue.prototype.getDirectedCoordinate=function(){return this._p1},Ue.prototype.getDx=function(){return this._dx},Ue.prototype.getLabel=function(){return this._label},Ue.prototype.getEdge=function(){return this._edge},Ue.prototype.getQuadrant=function(){return this._quadrant},Ue.prototype.getNode=function(){return this._node},Ue.prototype.toString=function(){var t=Math.atan2(this._dy,this._dx),e=this.getClass().getName(),n=e.lastIndexOf(".");return"  "+e.substring(n+1)+": "+this._p0+" - "+this._p1+" "+this._quadrant+":"+t+"   "+this._label},Ue.prototype.computeLabel=function(t){},Ue.prototype.init=function(t,e){this._p0=t,this._p1=e,this._dx=e.x-t.x,this._dy=e.y-t.y,this._quadrant=Be.quadrant(this._dx,this._dy),et.isTrue(!(0===this._dx&&0===this._dy),"EdgeEnd with identical endpoints found")},Ue.prototype.interfaces_=function(){return[E]},Ue.prototype.getClass=function(){return Ue};var ze=function(t){function e(){var e=arguments[0],n=arguments[1];if(t.call(this,e),this._isForward=null,this._isInResult=!1,this._isVisited=!1,this._sym=null,this._next=null,this._nextMin=null,this._edgeRing=null,this._minEdgeRing=null,this._depth=[0,-999,-999],this._isForward=n,n)this.init(e.getCoordinate(0),e.getCoordinate(1));else{var i=e.getNumPoints()-1;this.init(e.getCoordinate(i),e.getCoordinate(i-1))}this.computeDirectedLabel()}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.getNextMin=function(){return this._nextMin},e.prototype.getDepth=function(t){return this._depth[t]},e.prototype.setVisited=function(t){this._isVisited=t},e.prototype.computeDirectedLabel=function(){this._label=new Pe(this._edge.getLabel()),this._isForward||this._label.flip()},e.prototype.getNext=function(){return this._next},e.prototype.setDepth=function(t,e){if(-999!==this._depth[t]&&this._depth[t]!==e)throw new we("assigned depths do not match",this.getCoordinate());this._depth[t]=e},e.prototype.isInteriorAreaEdge=function(){for(var t=!0,e=0;e<2;e++)this._label.isArea(e)&&this._label.getLocation(e,Se.LEFT)===w.INTERIOR&&this._label.getLocation(e,Se.RIGHT)===w.INTERIOR||(t=!1);return t},e.prototype.setNextMin=function(t){this._nextMin=t},e.prototype.print=function(e){t.prototype.print.call(this,e),e.print(" "+this._depth[Se.LEFT]+"/"+this._depth[Se.RIGHT]),e.print(" ("+this.getDepthDelta()+")"),this._isInResult&&e.print(" inResult")},e.prototype.setMinEdgeRing=function(t){this._minEdgeRing=t},e.prototype.isLineEdge=function(){var t=this._label.isLine(0)||this._label.isLine(1),e=!this._label.isArea(0)||this._label.allPositionsEqual(0,w.EXTERIOR),n=!this._label.isArea(1)||this._label.allPositionsEqual(1,w.EXTERIOR);return t&&e&&n},e.prototype.setEdgeRing=function(t){this._edgeRing=t},e.prototype.getMinEdgeRing=function(){return this._minEdgeRing},e.prototype.getDepthDelta=function(){var t=this._edge.getDepthDelta();return this._isForward||(t=-t),t},e.prototype.setInResult=function(t){this._isInResult=t},e.prototype.getSym=function(){return this._sym},e.prototype.isForward=function(){return this._isForward},e.prototype.getEdge=function(){return this._edge},e.prototype.printEdge=function(t){this.print(t),t.print(" "),this._isForward?this._edge.print(t):this._edge.printReverse(t)},e.prototype.setSym=function(t){this._sym=t},e.prototype.setVisitedEdge=function(t){this.setVisited(t),this._sym.setVisited(t)},e.prototype.setEdgeDepths=function(t,e){var n=this.getEdge().getDepthDelta();this._isForward||(n=-n);var i=1;t===Se.LEFT&&(i=-1);var r=Se.opposite(t),o=e+n*i;this.setDepth(t,e),this.setDepth(r,o)},e.prototype.getEdgeRing=function(){return this._edgeRing},e.prototype.isInResult=function(){return this._isInResult},e.prototype.setNext=function(t){this._next=t},e.prototype.isVisited=function(){return this._isVisited},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e.depthFactor=function(t,e){return t===w.EXTERIOR&&e===w.INTERIOR?1:t===w.INTERIOR&&e===w.EXTERIOR?-1:0},e}(Ue),Xe=function(){};Xe.prototype.createNode=function(t){return new Ge(t,null)},Xe.prototype.interfaces_=function(){return[]},Xe.prototype.getClass=function(){return Xe};var Ye=function(){if(this._edges=new Nt,this._nodes=null,this._edgeEndList=new Nt,0===arguments.length)this._nodes=new qe(new Xe);else if(1===arguments.length){var t=arguments[0];this._nodes=new qe(t)}};Ye.prototype.printEdges=function(t){t.println("Edges:");for(var e=0;e<this._edges.size();e++){t.println("edge "+e+":");var n=this._edges.get(e);n.print(t),n.eiList.print(t)}},Ye.prototype.find=function(t){return this._nodes.find(t)},Ye.prototype.addNode=function(){if(arguments[0]instanceof Ge){var t=arguments[0];return this._nodes.addNode(t)}if(arguments[0]instanceof C){var e=arguments[0];return this._nodes.addNode(e)}},Ye.prototype.getNodeIterator=function(){return this._nodes.iterator()},Ye.prototype.linkResultDirectedEdges=function(){for(var t=this._nodes.iterator();t.hasNext();){t.next().getEdges().linkResultDirectedEdges()}},Ye.prototype.debugPrintln=function(t){Y.out.println(t)},Ye.prototype.isBoundaryNode=function(t,e){var n=this._nodes.find(e);if(null===n)return!1;var i=n.getLabel();return null!==i&&i.getLocation(t)===w.BOUNDARY},Ye.prototype.linkAllDirectedEdges=function(){for(var t=this._nodes.iterator();t.hasNext();){t.next().getEdges().linkAllDirectedEdges()}},Ye.prototype.matchInSameDirection=function(t,e,n,i){return!!t.equals(n)&&(at.computeOrientation(t,e,i)===at.COLLINEAR&&Be.quadrant(t,e)===Be.quadrant(n,i))},Ye.prototype.getEdgeEnds=function(){return this._edgeEndList},Ye.prototype.debugPrint=function(t){Y.out.print(t)},Ye.prototype.getEdgeIterator=function(){return this._edges.iterator()},Ye.prototype.findEdgeInSameDirection=function(t,e){for(var n=0;n<this._edges.size();n++){var i=this._edges.get(n),r=i.getCoordinates();if(this.matchInSameDirection(t,e,r[0],r[1]))return i;if(this.matchInSameDirection(t,e,r[r.length-1],r[r.length-2]))return i}return null},Ye.prototype.insertEdge=function(t){this._edges.add(t)},Ye.prototype.findEdgeEnd=function(t){for(var e=this.getEdgeEnds().iterator();e.hasNext();){var n=e.next();if(n.getEdge()===t)return n}return null},Ye.prototype.addEdges=function(t){for(var e=t.iterator();e.hasNext();){var n=e.next();this._edges.add(n);var i=new ze(n,!0),r=new ze(n,!1);i.setSym(r),r.setSym(i),this.add(i),this.add(r)}},Ye.prototype.add=function(t){this._nodes.add(t),this._edgeEndList.add(t)},Ye.prototype.getNodes=function(){return this._nodes.values()},Ye.prototype.findEdge=function(t,e){for(var n=0;n<this._edges.size();n++){var i=this._edges.get(n),r=i.getCoordinates();if(t.equals(r[0])&&e.equals(r[1]))return i}return null},Ye.prototype.interfaces_=function(){return[]},Ye.prototype.getClass=function(){return Ye},Ye.linkResultDirectedEdges=function(t){for(var e=t.iterator();e.hasNext();){e.next().getEdges().linkResultDirectedEdges()}};var ke=function(){this._geometryFactory=null,this._shellList=new Nt;var t=arguments[0];this._geometryFactory=t};ke.prototype.sortShellsAndHoles=function(t,e,n){for(var i=t.iterator();i.hasNext();){var r=i.next();r.isHole()?n.add(r):e.add(r)}},ke.prototype.computePolygons=function(t){for(var e=new Nt,n=t.iterator();n.hasNext();){var i=n.next().toPolygon(this._geometryFactory);e.add(i)}return e},ke.prototype.placeFreeHoles=function(t,e){for(var n=e.iterator();n.hasNext();){var i=n.next();if(null===i.getShell()){var r=this.findEdgeRingContaining(i,t);if(null===r)throw new we("unable to assign hole to a shell",i.getCoordinate(0));i.setShell(r)}}},ke.prototype.buildMinimalEdgeRings=function(t,e,n){for(var i=new Nt,r=t.iterator();r.hasNext();){var o=r.next();if(o.getMaxNodeDegree()>2){o.linkDirectedEdgesForMinimalEdgeRings();var s=o.buildMinimalRings(),a=this.findShell(s);null!==a?(this.placePolygonHoles(a,s),e.add(a)):n.addAll(s)}else i.add(o)}return i},ke.prototype.containsPoint=function(t){for(var e=this._shellList.iterator();e.hasNext();){if(e.next().containsPoint(t))return!0}return!1},ke.prototype.buildMaximalEdgeRings=function(t){for(var e=new Nt,n=t.iterator();n.hasNext();){var i=n.next();if(i.isInResult()&&i.getLabel().isArea()&&null===i.getEdgeRing()){var r=new Ae(i,this._geometryFactory);e.add(r),r.setInResult()}}return e},ke.prototype.placePolygonHoles=function(t,e){for(var n=e.iterator();n.hasNext();){var i=n.next();i.isHole()&&i.setShell(t)}},ke.prototype.getPolygons=function(){return this.computePolygons(this._shellList)},ke.prototype.findEdgeRingContaining=function(t,e){for(var n=t.getLinearRing(),i=n.getEnvelopeInternal(),r=n.getCoordinateN(0),o=null,s=null,a=e.iterator();a.hasNext();){var u=a.next(),l=u.getLinearRing(),c=l.getEnvelopeInternal();null!==o&&(s=o.getLinearRing().getEnvelopeInternal());var p=!1;c.contains(i)&&at.isPointInRing(r,l.getCoordinates())&&(p=!0),p&&(null===o||s.contains(c))&&(o=u)}return o},ke.prototype.findShell=function(t){for(var e=0,n=null,i=t.iterator();i.hasNext();){var r=i.next();r.isHole()||(n=r,e++)}return et.isTrue(e<=1,"found two shells in MinimalEdgeRing list"),n},ke.prototype.add=function(){if(1===arguments.length){var t=arguments[0];this.add(t.getEdgeEnds(),t.getNodes())}else if(2===arguments.length){var e=arguments[0],n=arguments[1];Ye.linkResultDirectedEdges(n);var i=this.buildMaximalEdgeRings(e),r=new Nt,o=this.buildMinimalEdgeRings(i,this._shellList,r);this.sortShellsAndHoles(o,this._shellList,r),this.placeFreeHoles(this._shellList,r)}},ke.prototype.interfaces_=function(){return[]},ke.prototype.getClass=function(){return ke};var je=function(){};je.prototype.getBounds=function(){},je.prototype.interfaces_=function(){return[]},je.prototype.getClass=function(){return je};var He=function(){this._bounds=null,this._item=null;var t=arguments[0],e=arguments[1];this._bounds=t,this._item=e};He.prototype.getItem=function(){return this._item},He.prototype.getBounds=function(){return this._bounds},He.prototype.interfaces_=function(){return[je,e]},He.prototype.getClass=function(){return He};var We=function(){this._size=null,this._items=null,this._size=0,this._items=new Nt,this._items.add(null)};We.prototype.poll=function(){if(this.isEmpty())return null;var t=this._items.get(1);return this._items.set(1,this._items.get(this._size)),this._size-=1,this.reorder(1),t},We.prototype.size=function(){return this._size},We.prototype.reorder=function(t){for(var e=null,n=this._items.get(t);2*t<=this._size&&((e=2*t)!==this._size&&this._items.get(e+1).compareTo(this._items.get(e))<0&&e++,this._items.get(e).compareTo(n)<0);t=e)this._items.set(t,this._items.get(e));this._items.set(t,n)},We.prototype.clear=function(){this._size=0,this._items.clear()},We.prototype.isEmpty=function(){return 0===this._size},We.prototype.add=function(t){this._items.add(null),this._size+=1;var e=this._size;for(this._items.set(0,t);t.compareTo(this._items.get(Math.trunc(e/2)))<0;e/=2)this._items.set(e,this._items.get(Math.trunc(e/2)));this._items.set(e,t)},We.prototype.interfaces_=function(){return[]},We.prototype.getClass=function(){return We};var Ke=function(){};Ke.prototype.visitItem=function(t){},Ke.prototype.interfaces_=function(){return[]},Ke.prototype.getClass=function(){return Ke};var Je=function(){};Je.prototype.insert=function(t,e){},Je.prototype.remove=function(t,e){},Je.prototype.query=function(){},Je.prototype.interfaces_=function(){return[]},Je.prototype.getClass=function(){return Je};var Qe=function(){if(this._childBoundables=new Nt,this._bounds=null,this._level=null,0===arguments.length);else if(1===arguments.length){var t=arguments[0];this._level=t}},Ze={serialVersionUID:{configurable:!0}};Qe.prototype.getLevel=function(){return this._level},Qe.prototype.size=function(){return this._childBoundables.size()},Qe.prototype.getChildBoundables=function(){return this._childBoundables},Qe.prototype.addChildBoundable=function(t){et.isTrue(null===this._bounds),this._childBoundables.add(t)},Qe.prototype.isEmpty=function(){return this._childBoundables.isEmpty()},Qe.prototype.getBounds=function(){return null===this._bounds&&(this._bounds=this.computeBounds()),this._bounds},Qe.prototype.interfaces_=function(){return[je,e]},Qe.prototype.getClass=function(){return Qe},Ze.serialVersionUID.get=function(){return 0x5a1e55ec41369800},Object.defineProperties(Qe,Ze);var $e=function(){};$e.reverseOrder=function(){return{compare:function(t,e){return e.compareTo(t)}}},$e.min=function(t){return $e.sort(t),t.get(0)},$e.sort=function(t,e){var n=t.toArray();e?Gt.sort(n,e):Gt.sort(n);for(var i=t.iterator(),r=0,o=n.length;r<o;r++)i.next(),i.set(n[r])},$e.singletonList=function(t){var e=new Nt;return e.add(t),e};var tn=function(){this._boundable1=null,this._boundable2=null,this._distance=null,this._itemDistance=null;var t=arguments[0],e=arguments[1],n=arguments[2];this._boundable1=t,this._boundable2=e,this._itemDistance=n,this._distance=this.distance()};tn.prototype.expandToQueue=function(t,e){var n=tn.isComposite(this._boundable1),i=tn.isComposite(this._boundable2);if(n&&i)return tn.area(this._boundable1)>tn.area(this._boundable2)?(this.expand(this._boundable1,this._boundable2,t,e),null):(this.expand(this._boundable2,this._boundable1,t,e),null);if(n)return this.expand(this._boundable1,this._boundable2,t,e),null;if(i)return this.expand(this._boundable2,this._boundable1,t,e),null;throw new m("neither boundable is composite")},tn.prototype.isLeaves=function(){return!(tn.isComposite(this._boundable1)||tn.isComposite(this._boundable2))},tn.prototype.compareTo=function(t){var e=t;return this._distance<e._distance?-1:this._distance>e._distance?1:0},tn.prototype.expand=function(t,e,n,i){for(var r=t.getChildBoundables().iterator();r.hasNext();){var o=r.next(),s=new tn(o,e,this._itemDistance);s.getDistance()<i&&n.add(s)}},tn.prototype.getBoundable=function(t){return 0===t?this._boundable1:this._boundable2},tn.prototype.getDistance=function(){return this._distance},tn.prototype.distance=function(){return this.isLeaves()?this._itemDistance.distance(this._boundable1,this._boundable2):this._boundable1.getBounds().distance(this._boundable2.getBounds())},tn.prototype.interfaces_=function(){return[E]},tn.prototype.getClass=function(){return tn},tn.area=function(t){return t.getBounds().getArea()},tn.isComposite=function(t){return t instanceof Qe};var en=function t(){if(this._root=null,this._built=!1,this._itemBoundables=new Nt,this._nodeCapacity=null,0===arguments.length){var e=t.DEFAULT_NODE_CAPACITY;this._nodeCapacity=e}else if(1===arguments.length){var n=arguments[0];et.isTrue(n>1,"Node capacity must be greater than 1"),this._nodeCapacity=n}},nn={IntersectsOp:{configurable:!0},serialVersionUID:{configurable:!0},DEFAULT_NODE_CAPACITY:{configurable:!0}};en.prototype.getNodeCapacity=function(){return this._nodeCapacity},en.prototype.lastNode=function(t){return t.get(t.size()-1)},en.prototype.size=function(){if(0===arguments.length)return this.isEmpty()?0:(this.build(),this.size(this._root));if(1===arguments.length){for(var t=0,e=arguments[0].getChildBoundables().iterator();e.hasNext();){var n=e.next();n instanceof Qe?t+=this.size(n):n instanceof He&&(t+=1)}return t}},en.prototype.removeItem=function(t,e){for(var n=null,i=t.getChildBoundables().iterator();i.hasNext();){var r=i.next();r instanceof He&&r.getItem()===e&&(n=r)}return null!==n&&(t.getChildBoundables().remove(n),!0)},en.prototype.itemsTree=function(){if(0===arguments.length){this.build();var t=this.itemsTree(this._root);return null===t?new Nt:t}if(1===arguments.length){for(var e=arguments[0],n=new Nt,i=e.getChildBoundables().iterator();i.hasNext();){var r=i.next();if(r instanceof Qe){var o=this.itemsTree(r);null!==o&&n.add(o)}else r instanceof He?n.add(r.getItem()):et.shouldNeverReachHere()}return n.size()<=0?null:n}},en.prototype.insert=function(t,e){et.isTrue(!this._built,"Cannot insert items into an STR packed R-tree after it has been built."),this._itemBoundables.add(new He(t,e))},en.prototype.boundablesAtLevel=function(){if(1===arguments.length){var t=arguments[0],e=new Nt;return this.boundablesAtLevel(t,this._root,e),e}if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2];if(et.isTrue(n>-2),i.getLevel()===n)return r.add(i),null;for(var o=i.getChildBoundables().iterator();o.hasNext();){var s=o.next();s instanceof Qe?this.boundablesAtLevel(n,s,r):(et.isTrue(s instanceof He),-1===n&&r.add(s))}return null}},en.prototype.query=function(){if(1===arguments.length){var t=arguments[0];this.build();var e=new Nt;return this.isEmpty()?e:(this.getIntersectsOp().intersects(this._root.getBounds(),t)&&this.query(t,this._root,e),e)}if(2===arguments.length){var n=arguments[0],i=arguments[1];if(this.build(),this.isEmpty())return null;this.getIntersectsOp().intersects(this._root.getBounds(),n)&&this.query(n,this._root,i)}else if(3===arguments.length)if(T(arguments[2],Ke)&&arguments[0]instanceof Object&&arguments[1]instanceof Qe)for(var r=arguments[0],o=arguments[1],s=arguments[2],a=o.getChildBoundables(),u=0;u<a.size();u++){var l=a.get(u);this.getIntersectsOp().intersects(l.getBounds(),r)&&(l instanceof Qe?this.query(r,l,s):l instanceof He?s.visitItem(l.getItem()):et.shouldNeverReachHere())}else if(T(arguments[2],xt)&&arguments[0]instanceof Object&&arguments[1]instanceof Qe)for(var c=arguments[0],p=arguments[1],h=arguments[2],f=p.getChildBoundables(),g=0;g<f.size();g++){var d=f.get(g);this.getIntersectsOp().intersects(d.getBounds(),c)&&(d instanceof Qe?this.query(c,d,h):d instanceof He?h.add(d.getItem()):et.shouldNeverReachHere())}},en.prototype.build=function(){if(this._built)return null;this._root=this._itemBoundables.isEmpty()?this.createNode(0):this.createHigherLevels(this._itemBoundables,-1),this._itemBoundables=null,this._built=!0},en.prototype.getRoot=function(){return this.build(),this._root},en.prototype.remove=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];return this.build(),!!this.getIntersectsOp().intersects(this._root.getBounds(),t)&&this.remove(t,this._root,e)}if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2],o=this.removeItem(i,r);if(o)return!0;for(var s=null,a=i.getChildBoundables().iterator();a.hasNext();){var u=a.next();if(this.getIntersectsOp().intersects(u.getBounds(),n)&&(u instanceof Qe&&(o=this.remove(n,u,r)))){s=u;break}}return null!==s&&s.getChildBoundables().isEmpty()&&i.getChildBoundables().remove(s),o}},en.prototype.createHigherLevels=function(t,e){et.isTrue(!t.isEmpty());var n=this.createParentBoundables(t,e+1);return 1===n.size()?n.get(0):this.createHigherLevels(n,e+1)},en.prototype.depth=function(){if(0===arguments.length)return this.isEmpty()?0:(this.build(),this.depth(this._root));if(1===arguments.length){for(var t=0,e=arguments[0].getChildBoundables().iterator();e.hasNext();){var n=e.next();if(n instanceof Qe){var i=this.depth(n);i>t&&(t=i)}}return t+1}},en.prototype.createParentBoundables=function(t,e){et.isTrue(!t.isEmpty());var n=new Nt;n.add(this.createNode(e));var i=new Nt(t);$e.sort(i,this.getComparator());for(var r=i.iterator();r.hasNext();){var o=r.next();this.lastNode(n).getChildBoundables().size()===this.getNodeCapacity()&&n.add(this.createNode(e)),this.lastNode(n).addChildBoundable(o)}return n},en.prototype.isEmpty=function(){return this._built?this._root.isEmpty():this._itemBoundables.isEmpty()},en.prototype.interfaces_=function(){return[e]},en.prototype.getClass=function(){return en},en.compareDoubles=function(t,e){return t>e?1:t<e?-1:0},nn.IntersectsOp.get=function(){return rn},nn.serialVersionUID.get=function(){return-0x35ef64c82d4c5400},nn.DEFAULT_NODE_CAPACITY.get=function(){return 10},Object.defineProperties(en,nn);var rn=function(){},on=function(){};on.prototype.distance=function(t,e){},on.prototype.interfaces_=function(){return[]},on.prototype.getClass=function(){return on};var sn=function(t){function n(e){e=e||n.DEFAULT_NODE_CAPACITY,t.call(this,e)}t&&(n.__proto__=t),(n.prototype=Object.create(t&&t.prototype)).constructor=n;var i={STRtreeNode:{configurable:!0},serialVersionUID:{configurable:!0},xComparator:{configurable:!0},yComparator:{configurable:!0},intersectsOp:{configurable:!0},DEFAULT_NODE_CAPACITY:{configurable:!0}};return n.prototype.createParentBoundablesFromVerticalSlices=function(t,e){et.isTrue(t.length>0);for(var n=new Nt,i=0;i<t.length;i++)n.addAll(this.createParentBoundablesFromVerticalSlice(t[i],e));return n},n.prototype.createNode=function(t){return new an(t)},n.prototype.size=function(){return 0===arguments.length?t.prototype.size.call(this):t.prototype.size.apply(this,arguments)},n.prototype.insert=function(){if(2!==arguments.length)return t.prototype.insert.apply(this,arguments);var e=arguments[0],n=arguments[1];if(e.isNull())return null;t.prototype.insert.call(this,e,n)},n.prototype.getIntersectsOp=function(){return n.intersectsOp},n.prototype.verticalSlices=function(t,e){for(var n=Math.trunc(Math.ceil(t.size()/e)),i=new Array(e).fill(null),r=t.iterator(),o=0;o<e;o++){i[o]=new Nt;for(var s=0;r.hasNext()&&s<n;){var a=r.next();i[o].add(a),s++}}return i},n.prototype.query=function(){if(1===arguments.length){var e=arguments[0];return t.prototype.query.call(this,e)}if(2===arguments.length){var n=arguments[0],i=arguments[1];t.prototype.query.call(this,n,i)}else if(3===arguments.length)if(T(arguments[2],Ke)&&arguments[0]instanceof Object&&arguments[1]instanceof Qe){var r=arguments[0],o=arguments[1],s=arguments[2];t.prototype.query.call(this,r,o,s)}else if(T(arguments[2],xt)&&arguments[0]instanceof Object&&arguments[1]instanceof Qe){var a=arguments[0],u=arguments[1],l=arguments[2];t.prototype.query.call(this,a,u,l)}},n.prototype.getComparator=function(){return n.yComparator},n.prototype.createParentBoundablesFromVerticalSlice=function(e,n){return t.prototype.createParentBoundables.call(this,e,n)},n.prototype.remove=function(){if(2===arguments.length){var e=arguments[0],n=arguments[1];return t.prototype.remove.call(this,e,n)}return t.prototype.remove.apply(this,arguments)},n.prototype.depth=function(){return 0===arguments.length?t.prototype.depth.call(this):t.prototype.depth.apply(this,arguments)},n.prototype.createParentBoundables=function(t,e){et.isTrue(!t.isEmpty());var i=Math.trunc(Math.ceil(t.size()/this.getNodeCapacity())),r=new Nt(t);$e.sort(r,n.xComparator);var o=this.verticalSlices(r,Math.trunc(Math.ceil(Math.sqrt(i))));return this.createParentBoundablesFromVerticalSlices(o,e)},n.prototype.nearestNeighbour=function(){if(1===arguments.length){if(T(arguments[0],on)){var t=arguments[0],e=new tn(this.getRoot(),this.getRoot(),t);return this.nearestNeighbour(e)}if(arguments[0]instanceof tn){var i=arguments[0];return this.nearestNeighbour(i,v.POSITIVE_INFINITY)}}else if(2===arguments.length){if(arguments[0]instanceof n&&T(arguments[1],on)){var r=arguments[0],o=arguments[1],s=new tn(this.getRoot(),r.getRoot(),o);return this.nearestNeighbour(s)}if(arguments[0]instanceof tn&&"number"==typeof arguments[1]){var a=arguments[0],u=arguments[1],l=null,c=new We;for(c.add(a);!c.isEmpty()&&u>0;){var p=c.poll(),h=p.getDistance();if(h>=u)break;p.isLeaves()?(u=h,l=p):p.expandToQueue(c,u)}return[l.getBoundable(0).getItem(),l.getBoundable(1).getItem()]}}else if(3===arguments.length){var f=arguments[0],g=arguments[1],d=arguments[2],y=new He(f,g),_=new tn(this.getRoot(),y,d);return this.nearestNeighbour(_)[0]}},n.prototype.interfaces_=function(){return[Je,e]},n.prototype.getClass=function(){return n},n.centreX=function(t){return n.avg(t.getMinX(),t.getMaxX())},n.avg=function(t,e){return(t+e)/2},n.centreY=function(t){return n.avg(t.getMinY(),t.getMaxY())},i.STRtreeNode.get=function(){return an},i.serialVersionUID.get=function(){return 0x39920f7d5f261e0},i.xComparator.get=function(){return{interfaces_:function(){return[N]},compare:function(e,i){return t.compareDoubles(n.centreX(e.getBounds()),n.centreX(i.getBounds()))}}},i.yComparator.get=function(){return{interfaces_:function(){return[N]},compare:function(e,i){return t.compareDoubles(n.centreY(e.getBounds()),n.centreY(i.getBounds()))}}},i.intersectsOp.get=function(){return{interfaces_:function(){return[t.IntersectsOp]},intersects:function(t,e){return t.intersects(e)}}},i.DEFAULT_NODE_CAPACITY.get=function(){return 10},Object.defineProperties(n,i),n}(en),an=function(t){function e(){var e=arguments[0];t.call(this,e)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.computeBounds=function(){for(var t=null,e=this.getChildBoundables().iterator();e.hasNext();){var n=e.next();null===t?t=new j(n.getBounds()):t.expandToInclude(n.getBounds())}return t},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Qe),un=function(){};un.prototype.interfaces_=function(){return[]},un.prototype.getClass=function(){return un},un.relativeSign=function(t,e){return t<e?-1:t>e?1:0},un.compare=function(t,e,n){if(e.equals2D(n))return 0;var i=un.relativeSign(e.x,n.x),r=un.relativeSign(e.y,n.y);switch(t){case 0:return un.compareValue(i,r);case 1:return un.compareValue(r,i);case 2:return un.compareValue(r,-i);case 3:return un.compareValue(-i,r);case 4:return un.compareValue(-i,-r);case 5:return un.compareValue(-r,-i);case 6:return un.compareValue(-r,i);case 7:return un.compareValue(i,-r)}return et.shouldNeverReachHere("invalid octant value"),0},un.compareValue=function(t,e){return t<0?-1:t>0?1:e<0?-1:e>0?1:0};var ln=function(){this._segString=null,this.coord=null,this.segmentIndex=null,this._segmentOctant=null,this._isInterior=null;var t=arguments[0],e=arguments[1],n=arguments[2],i=arguments[3];this._segString=t,this.coord=new C(e),this.segmentIndex=n,this._segmentOctant=i,this._isInterior=!e.equals2D(t.getCoordinate(n))};ln.prototype.getCoordinate=function(){return this.coord},ln.prototype.print=function(t){t.print(this.coord),t.print(" seg # = "+this.segmentIndex)},ln.prototype.compareTo=function(t){var e=t;return this.segmentIndex<e.segmentIndex?-1:this.segmentIndex>e.segmentIndex?1:this.coord.equals2D(e.coord)?0:un.compare(this._segmentOctant,this.coord,e.coord)},ln.prototype.isEndPoint=function(t){return 0===this.segmentIndex&&!this._isInterior||this.segmentIndex===t},ln.prototype.isInterior=function(){return this._isInterior},ln.prototype.interfaces_=function(){return[E]},ln.prototype.getClass=function(){return ln};var cn=function(){this._nodeMap=new p,this._edge=null;var t=arguments[0];this._edge=t};cn.prototype.getSplitCoordinates=function(){var t=new St;this.addEndpoints();for(var e=this.iterator(),n=e.next();e.hasNext();){var i=e.next();this.addEdgeCoordinates(n,i,t),n=i}return t.toCoordinateArray()},cn.prototype.addCollapsedNodes=function(){var t=new Nt;this.findCollapsesFromInsertedNodes(t),this.findCollapsesFromExistingVertices(t);for(var e=t.iterator();e.hasNext();){var n=e.next().intValue();this.add(this._edge.getCoordinate(n),n)}},cn.prototype.print=function(t){t.println("Intersections:");for(var e=this.iterator();e.hasNext();){e.next().print(t)}},cn.prototype.findCollapsesFromExistingVertices=function(t){for(var e=0;e<this._edge.size()-2;e++){var n=this._edge.getCoordinate(e),i=this._edge.getCoordinate(e+2);n.equals2D(i)&&t.add(new M(e+1))}},cn.prototype.addEdgeCoordinates=function(t,e,n){var i=this._edge.getCoordinate(e.segmentIndex),r=e.isInterior()||!e.coord.equals2D(i);n.add(new C(t.coord),!1);for(var o=t.segmentIndex+1;o<=e.segmentIndex;o++)n.add(this._edge.getCoordinate(o));r&&n.add(new C(e.coord))},cn.prototype.iterator=function(){return this._nodeMap.values().iterator()},cn.prototype.addSplitEdges=function(t){this.addEndpoints(),this.addCollapsedNodes();for(var e=this.iterator(),n=e.next();e.hasNext();){var i=e.next(),r=this.createSplitEdge(n,i);t.add(r),n=i}},cn.prototype.findCollapseIndex=function(t,e,n){if(!t.coord.equals2D(e.coord))return!1;var i=e.segmentIndex-t.segmentIndex;return e.isInterior()||i--,1===i&&(n[0]=t.segmentIndex+1,!0)},cn.prototype.findCollapsesFromInsertedNodes=function(t){for(var e=new Array(1).fill(null),n=this.iterator(),i=n.next();n.hasNext();){var r=n.next();this.findCollapseIndex(i,r,e)&&t.add(new M(e[0])),i=r}},cn.prototype.getEdge=function(){return this._edge},cn.prototype.addEndpoints=function(){var t=this._edge.size()-1;this.add(this._edge.getCoordinate(0),0),this.add(this._edge.getCoordinate(t),t)},cn.prototype.createSplitEdge=function(t,e){var n=e.segmentIndex-t.segmentIndex+2,i=this._edge.getCoordinate(e.segmentIndex),r=e.isInterior()||!e.coord.equals2D(i);r||n--;var o=new Array(n).fill(null),s=0;o[s++]=new C(t.coord);for(var a=t.segmentIndex+1;a<=e.segmentIndex;a++)o[s++]=this._edge.getCoordinate(a);return r&&(o[s]=new C(e.coord)),new gn(o,this._edge.getData())},cn.prototype.add=function(t,e){var n=new ln(this._edge,t,e,this._edge.getSegmentOctant(e)),i=this._nodeMap.get(n);return null!==i?(et.isTrue(i.coord.equals2D(t),"Found equal nodes with different coordinates"),i):(this._nodeMap.put(n,n),n)},cn.prototype.checkSplitEdgesCorrectness=function(t){var e=this._edge.getCoordinates(),n=t.get(0).getCoordinate(0);if(!n.equals2D(e[0]))throw new $("bad split edge start point at "+n);var i=t.get(t.size()-1).getCoordinates(),r=i[i.length-1];if(!r.equals2D(e[e.length-1]))throw new $("bad split edge end point at "+r)},cn.prototype.interfaces_=function(){return[]},cn.prototype.getClass=function(){return cn};var pn=function(){};pn.prototype.interfaces_=function(){return[]},pn.prototype.getClass=function(){return pn},pn.octant=function(){if("number"==typeof arguments[0]&&"number"==typeof arguments[1]){var t=arguments[0],e=arguments[1];if(0===t&&0===e)throw new m("Cannot compute the octant for point ( "+t+", "+e+" )");var n=Math.abs(t),i=Math.abs(e);return t>=0?e>=0?n>=i?0:1:n>=i?7:6:e>=0?n>=i?3:2:n>=i?4:5}if(arguments[0]instanceof C&&arguments[1]instanceof C){var r=arguments[0],o=arguments[1],s=o.x-r.x,a=o.y-r.y;if(0===s&&0===a)throw new m("Cannot compute the octant for two identical points "+r);return pn.octant(s,a)}};var hn=function(){};hn.prototype.getCoordinates=function(){},hn.prototype.size=function(){},hn.prototype.getCoordinate=function(t){},hn.prototype.isClosed=function(){},hn.prototype.setData=function(t){},hn.prototype.getData=function(){},hn.prototype.interfaces_=function(){return[]},hn.prototype.getClass=function(){return hn};var fn=function(){};fn.prototype.addIntersection=function(t,e){},fn.prototype.interfaces_=function(){return[hn]},fn.prototype.getClass=function(){return fn};var gn=function(){this._nodeList=new cn(this),this._pts=null,this._data=null;var t=arguments[0],e=arguments[1];this._pts=t,this._data=e};gn.prototype.getCoordinates=function(){return this._pts},gn.prototype.size=function(){return this._pts.length},gn.prototype.getCoordinate=function(t){return this._pts[t]},gn.prototype.isClosed=function(){return this._pts[0].equals(this._pts[this._pts.length-1])},gn.prototype.getSegmentOctant=function(t){return t===this._pts.length-1?-1:this.safeOctant(this.getCoordinate(t),this.getCoordinate(t+1))},gn.prototype.setData=function(t){this._data=t},gn.prototype.safeOctant=function(t,e){return t.equals2D(e)?0:pn.octant(t,e)},gn.prototype.getData=function(){return this._data},gn.prototype.addIntersection=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];this.addIntersectionNode(t,e)}else if(4===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[3],o=new C(n.getIntersection(r));this.addIntersection(o,i)}},gn.prototype.toString=function(){return Z.toLineString(new ue(this._pts))},gn.prototype.getNodeList=function(){return this._nodeList},gn.prototype.addIntersectionNode=function(t,e){var n=e,i=n+1;if(i<this._pts.length){var r=this._pts[i];t.equals2D(r)&&(n=i)}return this._nodeList.add(t,n)},gn.prototype.addIntersections=function(t,e,n){for(var i=0;i<t.getIntersectionNum();i++)this.addIntersection(t,e,n,i)},gn.prototype.interfaces_=function(){return[fn]},gn.prototype.getClass=function(){return gn},gn.getNodedSubstrings=function(){if(1===arguments.length){var t=arguments[0],e=new Nt;return gn.getNodedSubstrings(t,e),e}if(2===arguments.length)for(var n=arguments[0],i=arguments[1],r=n.iterator();r.hasNext();){r.next().getNodeList().addSplitEdges(i)}};var dn=function(){if(this.p0=null,this.p1=null,0===arguments.length)this.p0=new C,this.p1=new C;else if(1===arguments.length){var t=arguments[0];this.p0=new C(t.p0),this.p1=new C(t.p1)}else if(2===arguments.length)this.p0=arguments[0],this.p1=arguments[1];else if(4===arguments.length){var e=arguments[0],n=arguments[1],i=arguments[2],r=arguments[3];this.p0=new C(e,n),this.p1=new C(i,r)}},yn={serialVersionUID:{configurable:!0}};dn.prototype.minX=function(){return Math.min(this.p0.x,this.p1.x)},dn.prototype.orientationIndex=function(){if(arguments[0]instanceof dn){var t=arguments[0],e=at.orientationIndex(this.p0,this.p1,t.p0),n=at.orientationIndex(this.p0,this.p1,t.p1);return e>=0&&n>=0?Math.max(e,n):e<=0&&n<=0?Math.max(e,n):0}if(arguments[0]instanceof C){var i=arguments[0];return at.orientationIndex(this.p0,this.p1,i)}},dn.prototype.toGeometry=function(t){return t.createLineString([this.p0,this.p1])},dn.prototype.isVertical=function(){return this.p0.x===this.p1.x},dn.prototype.equals=function(t){if(!(t instanceof dn))return!1;var e=t;return this.p0.equals(e.p0)&&this.p1.equals(e.p1)},dn.prototype.intersection=function(t){var e=new rt;return e.computeIntersection(this.p0,this.p1,t.p0,t.p1),e.hasIntersection()?e.getIntersection(0):null},dn.prototype.project=function(){if(arguments[0]instanceof C){var t=arguments[0];if(t.equals(this.p0)||t.equals(this.p1))return new C(t);var e=this.projectionFactor(t),n=new C;return n.x=this.p0.x+e*(this.p1.x-this.p0.x),n.y=this.p0.y+e*(this.p1.y-this.p0.y),n}if(arguments[0]instanceof dn){var i=arguments[0],r=this.projectionFactor(i.p0),o=this.projectionFactor(i.p1);if(r>=1&&o>=1)return null;if(r<=0&&o<=0)return null;var s=this.project(i.p0);r<0&&(s=this.p0),r>1&&(s=this.p1);var a=this.project(i.p1);return o<0&&(a=this.p0),o>1&&(a=this.p1),new dn(s,a)}},dn.prototype.normalize=function(){this.p1.compareTo(this.p0)<0&&this.reverse()},dn.prototype.angle=function(){return Math.atan2(this.p1.y-this.p0.y,this.p1.x-this.p0.x)},dn.prototype.getCoordinate=function(t){return 0===t?this.p0:this.p1},dn.prototype.distancePerpendicular=function(t){return at.distancePointLinePerpendicular(t,this.p0,this.p1)},dn.prototype.minY=function(){return Math.min(this.p0.y,this.p1.y)},dn.prototype.midPoint=function(){return dn.midPoint(this.p0,this.p1)},dn.prototype.projectionFactor=function(t){if(t.equals(this.p0))return 0;if(t.equals(this.p1))return 1;var e=this.p1.x-this.p0.x,n=this.p1.y-this.p0.y,i=e*e+n*n;if(i<=0)return v.NaN;return((t.x-this.p0.x)*e+(t.y-this.p0.y)*n)/i},dn.prototype.closestPoints=function(t){var e=this.intersection(t);if(null!==e)return[e,e];var n=new Array(2).fill(null),i=v.MAX_VALUE,r=null,o=this.closestPoint(t.p0);i=o.distance(t.p0),n[0]=o,n[1]=t.p0;var s=this.closestPoint(t.p1);(r=s.distance(t.p1))<i&&(i=r,n[0]=s,n[1]=t.p1);var a=t.closestPoint(this.p0);(r=a.distance(this.p0))<i&&(i=r,n[0]=this.p0,n[1]=a);var u=t.closestPoint(this.p1);return(r=u.distance(this.p1))<i&&(i=r,n[0]=this.p1,n[1]=u),n},dn.prototype.closestPoint=function(t){var e=this.projectionFactor(t);if(e>0&&e<1)return this.project(t);return this.p0.distance(t)<this.p1.distance(t)?this.p0:this.p1},dn.prototype.maxX=function(){return Math.max(this.p0.x,this.p1.x)},dn.prototype.getLength=function(){return this.p0.distance(this.p1)},dn.prototype.compareTo=function(t){var e=t,n=this.p0.compareTo(e.p0);return 0!==n?n:this.p1.compareTo(e.p1)},dn.prototype.reverse=function(){var t=this.p0;this.p0=this.p1,this.p1=t},dn.prototype.equalsTopo=function(t){return this.p0.equals(t.p0)&&(this.p1.equals(t.p1)||this.p0.equals(t.p1))&&this.p1.equals(t.p0)},dn.prototype.lineIntersection=function(t){try{return k.intersection(this.p0,this.p1,t.p0,t.p1)}catch(t){if(!(t instanceof X))throw t}return null},dn.prototype.maxY=function(){return Math.max(this.p0.y,this.p1.y)},dn.prototype.pointAlongOffset=function(t,e){var n=this.p0.x+t*(this.p1.x-this.p0.x),i=this.p0.y+t*(this.p1.y-this.p0.y),r=this.p1.x-this.p0.x,o=this.p1.y-this.p0.y,s=Math.sqrt(r*r+o*o),a=0,u=0;if(0!==e){if(s<=0)throw new Error("Cannot compute offset from zero-length line segment");a=e*r/s,u=e*o/s}return new C(n-u,i+a)},dn.prototype.setCoordinates=function(){if(1===arguments.length){var t=arguments[0];this.setCoordinates(t.p0,t.p1)}else if(2===arguments.length){var e=arguments[0],n=arguments[1];this.p0.x=e.x,this.p0.y=e.y,this.p1.x=n.x,this.p1.y=n.y}},dn.prototype.segmentFraction=function(t){var e=this.projectionFactor(t);return e<0?e=0:(e>1||v.isNaN(e))&&(e=1),e},dn.prototype.toString=function(){return"LINESTRING( "+this.p0.x+" "+this.p0.y+", "+this.p1.x+" "+this.p1.y+")"},dn.prototype.isHorizontal=function(){return this.p0.y===this.p1.y},dn.prototype.distance=function(){if(arguments[0]instanceof dn){var t=arguments[0];return at.distanceLineLine(this.p0,this.p1,t.p0,t.p1)}if(arguments[0]instanceof C){var e=arguments[0];return at.distancePointLine(e,this.p0,this.p1)}},dn.prototype.pointAlong=function(t){var e=new C;return e.x=this.p0.x+t*(this.p1.x-this.p0.x),e.y=this.p0.y+t*(this.p1.y-this.p0.y),e},dn.prototype.hashCode=function(){var t=v.doubleToLongBits(this.p0.x);t^=31*v.doubleToLongBits(this.p0.y);var e=Math.trunc(t)^Math.trunc(t>>32),n=v.doubleToLongBits(this.p1.x);n^=31*v.doubleToLongBits(this.p1.y);return e^(Math.trunc(n)^Math.trunc(n>>32))},dn.prototype.interfaces_=function(){return[E,e]},dn.prototype.getClass=function(){return dn},dn.midPoint=function(t,e){return new C((t.x+e.x)/2,(t.y+e.y)/2)},yn.serialVersionUID.get=function(){return 0x2d2172135f411c00},Object.defineProperties(dn,yn);var _n=function(){this.tempEnv1=new j,this.tempEnv2=new j,this._overlapSeg1=new dn,this._overlapSeg2=new dn};_n.prototype.overlap=function(){if(2===arguments.length);else if(4===arguments.length){var t=arguments[0],e=arguments[1],n=arguments[2],i=arguments[3];t.getLineSegment(e,this._overlapSeg1),n.getLineSegment(i,this._overlapSeg2),this.overlap(this._overlapSeg1,this._overlapSeg2)}},_n.prototype.interfaces_=function(){return[]},_n.prototype.getClass=function(){return _n};var mn=function(){this._pts=null,this._start=null,this._end=null,this._env=null,this._context=null,this._id=null;var t=arguments[0],e=arguments[1],n=arguments[2],i=arguments[3];this._pts=t,this._start=e,this._end=n,this._context=i};mn.prototype.getLineSegment=function(t,e){e.p0=this._pts[t],e.p1=this._pts[t+1]},mn.prototype.computeSelect=function(t,e,n,i){var r=this._pts[e],o=this._pts[n];if(i.tempEnv1.init(r,o),n-e==1)return i.select(this,e),null;if(!t.intersects(i.tempEnv1))return null;var s=Math.trunc((e+n)/2);e<s&&this.computeSelect(t,e,s,i),s<n&&this.computeSelect(t,s,n,i)},mn.prototype.getCoordinates=function(){for(var t=new Array(this._end-this._start+1).fill(null),e=0,n=this._start;n<=this._end;n++)t[e++]=this._pts[n];return t},mn.prototype.computeOverlaps=function(t,e){this.computeOverlapsInternal(this._start,this._end,t,t._start,t._end,e)},mn.prototype.setId=function(t){this._id=t},mn.prototype.select=function(t,e){this.computeSelect(t,this._start,this._end,e)},mn.prototype.getEnvelope=function(){if(null===this._env){var t=this._pts[this._start],e=this._pts[this._end];this._env=new j(t,e)}return this._env},mn.prototype.getEndIndex=function(){return this._end},mn.prototype.getStartIndex=function(){return this._start},mn.prototype.getContext=function(){return this._context},mn.prototype.getId=function(){return this._id},mn.prototype.computeOverlapsInternal=function(t,e,n,i,r,o){var s=this._pts[t],a=this._pts[e],u=n._pts[i],l=n._pts[r];if(e-t==1&&r-i==1)return o.overlap(this,t,n,i),null;if(o.tempEnv1.init(s,a),o.tempEnv2.init(u,l),!o.tempEnv1.intersects(o.tempEnv2))return null;var c=Math.trunc((t+e)/2),p=Math.trunc((i+r)/2);t<c&&(i<p&&this.computeOverlapsInternal(t,c,n,i,p,o),p<r&&this.computeOverlapsInternal(t,c,n,p,r,o)),c<e&&(i<p&&this.computeOverlapsInternal(c,e,n,i,p,o),p<r&&this.computeOverlapsInternal(c,e,n,p,r,o))},mn.prototype.interfaces_=function(){return[]},mn.prototype.getClass=function(){return mn};var vn=function(){};vn.prototype.interfaces_=function(){return[]},vn.prototype.getClass=function(){return vn},vn.getChainStartIndices=function(t){var e=0,n=new Nt;n.add(new M(e));do{var i=vn.findChainEnd(t,e);n.add(new M(i)),e=i}while(e<t.length-1);return vn.toIntArray(n)},vn.findChainEnd=function(t,e){for(var n=e;n<t.length-1&&t[n].equals2D(t[n+1]);)n++;if(n>=t.length-1)return t.length-1;for(var i=Be.quadrant(t[n],t[n+1]),r=e+1;r<t.length;){if(!t[r-1].equals2D(t[r])){if(Be.quadrant(t[r-1],t[r])!==i)break}r++}return r-1},vn.getChains=function(){if(1===arguments.length){var t=arguments[0];return vn.getChains(t,null)}if(2===arguments.length){for(var e=arguments[0],n=arguments[1],i=new Nt,r=vn.getChainStartIndices(e),o=0;o<r.length-1;o++){var s=new mn(e,r[o],r[o+1],n);i.add(s)}return i}},vn.toIntArray=function(t){for(var e=new Array(t.size()).fill(null),n=0;n<e.length;n++)e[n]=t.get(n).intValue();return e};var In=function(){};In.prototype.computeNodes=function(t){},In.prototype.getNodedSubstrings=function(){},In.prototype.interfaces_=function(){return[]},In.prototype.getClass=function(){return In};var En=function(){if(this._segInt=null,0===arguments.length);else if(1===arguments.length){var t=arguments[0];this.setSegmentIntersector(t)}};En.prototype.setSegmentIntersector=function(t){this._segInt=t},En.prototype.interfaces_=function(){return[In]},En.prototype.getClass=function(){return En};var xn=function(t){function e(e){e?t.call(this,e):t.call(this),this._monoChains=new Nt,this._index=new sn,this._idCounter=0,this._nodedSegStrings=null,this._nOverlaps=0}t&&(e.__proto__=t),(e.prototype=Object.create(t&&t.prototype)).constructor=e;var n={SegmentOverlapAction:{configurable:!0}};return e.prototype.getMonotoneChains=function(){return this._monoChains},e.prototype.getNodedSubstrings=function(){return gn.getNodedSubstrings(this._nodedSegStrings)},e.prototype.getIndex=function(){return this._index},e.prototype.add=function(t){for(var e=vn.getChains(t.getCoordinates(),t).iterator();e.hasNext();){var n=e.next();n.setId(this._idCounter++),this._index.insert(n.getEnvelope(),n),this._monoChains.add(n)}},e.prototype.computeNodes=function(t){this._nodedSegStrings=t;for(var e=t.iterator();e.hasNext();)this.add(e.next());this.intersectChains()},e.prototype.intersectChains=function(){for(var t=new Nn(this._segInt),e=this._monoChains.iterator();e.hasNext();)for(var n=e.next(),i=this._index.query(n.getEnvelope()).iterator();i.hasNext();){var r=i.next();if(r.getId()>n.getId()&&(n.computeOverlaps(r,t),this._nOverlaps++),this._segInt.isDone())return null}},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},n.SegmentOverlapAction.get=function(){return Nn},Object.defineProperties(e,n),e}(En),Nn=function(t){function e(){t.call(this),this._si=null;var e=arguments[0];this._si=e}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.overlap=function(){if(4!==arguments.length)return t.prototype.overlap.apply(this,arguments);var e=arguments[0],n=arguments[1],i=arguments[2],r=arguments[3],o=e.getContext(),s=i.getContext();this._si.processIntersections(o,n,s,r)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(_n),Cn=function t(){if(this._quadrantSegments=t.DEFAULT_QUADRANT_SEGMENTS,this._endCapStyle=t.CAP_ROUND,this._joinStyle=t.JOIN_ROUND,this._mitreLimit=t.DEFAULT_MITRE_LIMIT,this._isSingleSided=!1,this._simplifyFactor=t.DEFAULT_SIMPLIFY_FACTOR,0===arguments.length);else if(1===arguments.length){var e=arguments[0];this.setQuadrantSegments(e)}else if(2===arguments.length){var n=arguments[0],i=arguments[1];this.setQuadrantSegments(n),this.setEndCapStyle(i)}else if(4===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2],a=arguments[3];this.setQuadrantSegments(r),this.setEndCapStyle(o),this.setJoinStyle(s),this.setMitreLimit(a)}},Sn={CAP_ROUND:{configurable:!0},CAP_FLAT:{configurable:!0},CAP_SQUARE:{configurable:!0},JOIN_ROUND:{configurable:!0},JOIN_MITRE:{configurable:!0},JOIN_BEVEL:{configurable:!0},DEFAULT_QUADRANT_SEGMENTS:{configurable:!0},DEFAULT_MITRE_LIMIT:{configurable:!0},DEFAULT_SIMPLIFY_FACTOR:{configurable:!0}};Cn.prototype.getEndCapStyle=function(){return this._endCapStyle},Cn.prototype.isSingleSided=function(){return this._isSingleSided},Cn.prototype.setQuadrantSegments=function(t){this._quadrantSegments=t,0===this._quadrantSegments&&(this._joinStyle=Cn.JOIN_BEVEL),this._quadrantSegments<0&&(this._joinStyle=Cn.JOIN_MITRE,this._mitreLimit=Math.abs(this._quadrantSegments)),t<=0&&(this._quadrantSegments=1),this._joinStyle!==Cn.JOIN_ROUND&&(this._quadrantSegments=Cn.DEFAULT_QUADRANT_SEGMENTS)},Cn.prototype.getJoinStyle=function(){return this._joinStyle},Cn.prototype.setJoinStyle=function(t){this._joinStyle=t},Cn.prototype.setSimplifyFactor=function(t){this._simplifyFactor=t<0?0:t},Cn.prototype.getSimplifyFactor=function(){return this._simplifyFactor},Cn.prototype.getQuadrantSegments=function(){return this._quadrantSegments},Cn.prototype.setEndCapStyle=function(t){this._endCapStyle=t},Cn.prototype.getMitreLimit=function(){return this._mitreLimit},Cn.prototype.setMitreLimit=function(t){this._mitreLimit=t},Cn.prototype.setSingleSided=function(t){this._isSingleSided=t},Cn.prototype.interfaces_=function(){return[]},Cn.prototype.getClass=function(){return Cn},Cn.bufferDistanceError=function(t){var e=Math.PI/2/t;return 1-Math.cos(e/2)},Sn.CAP_ROUND.get=function(){return 1},Sn.CAP_FLAT.get=function(){return 2},Sn.CAP_SQUARE.get=function(){return 3},Sn.JOIN_ROUND.get=function(){return 1},Sn.JOIN_MITRE.get=function(){return 2},Sn.JOIN_BEVEL.get=function(){return 3},Sn.DEFAULT_QUADRANT_SEGMENTS.get=function(){return 8},Sn.DEFAULT_MITRE_LIMIT.get=function(){return 5},Sn.DEFAULT_SIMPLIFY_FACTOR.get=function(){return.01},Object.defineProperties(Cn,Sn);var Ln=function(t){this._distanceTol=null,this._isDeleted=null,this._angleOrientation=at.COUNTERCLOCKWISE,this._inputLine=t||null},bn={INIT:{configurable:!0},DELETE:{configurable:!0},KEEP:{configurable:!0},NUM_PTS_TO_CHECK:{configurable:!0}};Ln.prototype.isDeletable=function(t,e,n,i){var r=this._inputLine[t],o=this._inputLine[e],s=this._inputLine[n];return!!this.isConcave(r,o,s)&&(!!this.isShallow(r,o,s,i)&&this.isShallowSampled(r,o,t,n,i))},Ln.prototype.deleteShallowConcavities=function(){for(var t=1,e=this.findNextNonDeletedIndex(t),n=this.findNextNonDeletedIndex(e),i=!1;n<this._inputLine.length;){var r=!1;this.isDeletable(t,e,n,this._distanceTol)&&(this._isDeleted[e]=Ln.DELETE,r=!0,i=!0),t=r?n:e,e=this.findNextNonDeletedIndex(t),n=this.findNextNonDeletedIndex(e)}return i},Ln.prototype.isShallowConcavity=function(t,e,n,i){if(!(at.computeOrientation(t,e,n)===this._angleOrientation))return!1;return at.distancePointLine(e,t,n)<i},Ln.prototype.isShallowSampled=function(t,e,n,i,r){var o=Math.trunc((i-n)/Ln.NUM_PTS_TO_CHECK);o<=0&&(o=1);for(var s=n;s<i;s+=o)if(!this.isShallow(t,e,this._inputLine[s],r))return!1;return!0},Ln.prototype.isConcave=function(t,e,n){var i=at.computeOrientation(t,e,n)===this._angleOrientation;return i},Ln.prototype.simplify=function(t){this._distanceTol=Math.abs(t),t<0&&(this._angleOrientation=at.CLOCKWISE),this._isDeleted=new Array(this._inputLine.length).fill(null);var e=!1;do{e=this.deleteShallowConcavities()}while(e);return this.collapseLine()},Ln.prototype.findNextNonDeletedIndex=function(t){for(var e=t+1;e<this._inputLine.length&&this._isDeleted[e]===Ln.DELETE;)e++;return e},Ln.prototype.isShallow=function(t,e,n,i){return at.distancePointLine(e,t,n)<i},Ln.prototype.collapseLine=function(){for(var t=new St,e=0;e<this._inputLine.length;e++)this._isDeleted[e]!==Ln.DELETE&&t.add(this._inputLine[e]);return t.toCoordinateArray()},Ln.prototype.interfaces_=function(){return[]},Ln.prototype.getClass=function(){return Ln},Ln.simplify=function(t,e){return new Ln(t).simplify(e)},bn.INIT.get=function(){return 0},bn.DELETE.get=function(){return 1},bn.KEEP.get=function(){return 1},bn.NUM_PTS_TO_CHECK.get=function(){return 10},Object.defineProperties(Ln,bn);var wn=function(){this._ptList=null,this._precisionModel=null,this._minimimVertexDistance=0,this._ptList=new Nt},On={COORDINATE_ARRAY_TYPE:{configurable:!0}};wn.prototype.getCoordinates=function(){return this._ptList.toArray(wn.COORDINATE_ARRAY_TYPE)},wn.prototype.setPrecisionModel=function(t){this._precisionModel=t},wn.prototype.addPt=function(t){var e=new C(t);if(this._precisionModel.makePrecise(e),this.isRedundant(e))return null;this._ptList.add(e)},wn.prototype.revere=function(){},wn.prototype.addPts=function(t,e){if(e)for(var n=0;n<t.length;n++)this.addPt(t[n]);else for(var i=t.length-1;i>=0;i--)this.addPt(t[i])},wn.prototype.isRedundant=function(t){if(this._ptList.size()<1)return!1;var e=this._ptList.get(this._ptList.size()-1);return t.distance(e)<this._minimimVertexDistance},wn.prototype.toString=function(){return(new _e).createLineString(this.getCoordinates()).toString()},wn.prototype.closeRing=function(){if(this._ptList.size()<1)return null;var t=new C(this._ptList.get(0)),e=this._ptList.get(this._ptList.size()-1);if(t.equals(e))return null;this._ptList.add(t)},wn.prototype.setMinimumVertexDistance=function(t){this._minimimVertexDistance=t},wn.prototype.interfaces_=function(){return[]},wn.prototype.getClass=function(){return wn},On.COORDINATE_ARRAY_TYPE.get=function(){return new Array(0).fill(null)},Object.defineProperties(wn,On);var Tn=function(){},Rn={PI_TIMES_2:{configurable:!0},PI_OVER_2:{configurable:!0},PI_OVER_4:{configurable:!0},COUNTERCLOCKWISE:{configurable:!0},CLOCKWISE:{configurable:!0},NONE:{configurable:!0}};Tn.prototype.interfaces_=function(){return[]},Tn.prototype.getClass=function(){return Tn},Tn.toDegrees=function(t){return 180*t/Math.PI},Tn.normalize=function(t){for(;t>Math.PI;)t-=Tn.PI_TIMES_2;for(;t<=-Math.PI;)t+=Tn.PI_TIMES_2;return t},Tn.angle=function(){if(1===arguments.length){var t=arguments[0];return Math.atan2(t.y,t.x)}if(2===arguments.length){var e=arguments[0],n=arguments[1],i=n.x-e.x,r=n.y-e.y;return Math.atan2(r,i)}},Tn.isAcute=function(t,e,n){var i=t.x-e.x,r=t.y-e.y;return i*(n.x-e.x)+r*(n.y-e.y)>0},Tn.isObtuse=function(t,e,n){var i=t.x-e.x,r=t.y-e.y;return i*(n.x-e.x)+r*(n.y-e.y)<0},Tn.interiorAngle=function(t,e,n){var i=Tn.angle(e,t),r=Tn.angle(e,n);return Math.abs(r-i)},Tn.normalizePositive=function(t){if(t<0){for(;t<0;)t+=Tn.PI_TIMES_2;t>=Tn.PI_TIMES_2&&(t=0)}else{for(;t>=Tn.PI_TIMES_2;)t-=Tn.PI_TIMES_2;t<0&&(t=0)}return t},Tn.angleBetween=function(t,e,n){var i=Tn.angle(e,t),r=Tn.angle(e,n);return Tn.diff(i,r)},Tn.diff=function(t,e){var n=null;return(n=t<e?e-t:t-e)>Math.PI&&(n=2*Math.PI-n),n},Tn.toRadians=function(t){return t*Math.PI/180},Tn.getTurn=function(t,e){var n=Math.sin(e-t);return n>0?Tn.COUNTERCLOCKWISE:n<0?Tn.CLOCKWISE:Tn.NONE},Tn.angleBetweenOriented=function(t,e,n){var i=Tn.angle(e,t),r=Tn.angle(e,n)-i;return r<=-Math.PI?r+Tn.PI_TIMES_2:r>Math.PI?r-Tn.PI_TIMES_2:r},Rn.PI_TIMES_2.get=function(){return 2*Math.PI},Rn.PI_OVER_2.get=function(){return Math.PI/2},Rn.PI_OVER_4.get=function(){return Math.PI/4},Rn.COUNTERCLOCKWISE.get=function(){return at.COUNTERCLOCKWISE},Rn.CLOCKWISE.get=function(){return at.CLOCKWISE},Rn.NONE.get=function(){return at.COLLINEAR},Object.defineProperties(Tn,Rn);var Pn=function t(){this._maxCurveSegmentError=0,this._filletAngleQuantum=null,this._closingSegLengthFactor=1,this._segList=null,this._distance=0,this._precisionModel=null,this._bufParams=null,this._li=null,this._s0=null,this._s1=null,this._s2=null,this._seg0=new dn,this._seg1=new dn,this._offset0=new dn,this._offset1=new dn,this._side=0,this._hasNarrowConcaveAngle=!1;var e=arguments[0],n=arguments[1],i=arguments[2];this._precisionModel=e,this._bufParams=n,this._li=new rt,this._filletAngleQuantum=Math.PI/2/n.getQuadrantSegments(),n.getQuadrantSegments()>=8&&n.getJoinStyle()===Cn.JOIN_ROUND&&(this._closingSegLengthFactor=t.MAX_CLOSING_SEG_LEN_FACTOR),this.init(i)},Dn={OFFSET_SEGMENT_SEPARATION_FACTOR:{configurable:!0},INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR:{configurable:!0},CURVE_VERTEX_SNAP_DISTANCE_FACTOR:{configurable:!0},MAX_CLOSING_SEG_LEN_FACTOR:{configurable:!0}};Pn.prototype.addNextSegment=function(t,e){if(this._s0=this._s1,this._s1=this._s2,this._s2=t,this._seg0.setCoordinates(this._s0,this._s1),this.computeOffsetSegment(this._seg0,this._side,this._distance,this._offset0),this._seg1.setCoordinates(this._s1,this._s2),this.computeOffsetSegment(this._seg1,this._side,this._distance,this._offset1),this._s1.equals(this._s2))return null;var n=at.computeOrientation(this._s0,this._s1,this._s2),i=n===at.CLOCKWISE&&this._side===Se.LEFT||n===at.COUNTERCLOCKWISE&&this._side===Se.RIGHT;0===n?this.addCollinear(e):i?this.addOutsideTurn(n,e):this.addInsideTurn(n,e)},Pn.prototype.addLineEndCap=function(t,e){var n=new dn(t,e),i=new dn;this.computeOffsetSegment(n,Se.LEFT,this._distance,i);var r=new dn;this.computeOffsetSegment(n,Se.RIGHT,this._distance,r);var o=e.x-t.x,s=e.y-t.y,a=Math.atan2(s,o);switch(this._bufParams.getEndCapStyle()){case Cn.CAP_ROUND:this._segList.addPt(i.p1),this.addFilletArc(e,a+Math.PI/2,a-Math.PI/2,at.CLOCKWISE,this._distance),this._segList.addPt(r.p1);break;case Cn.CAP_FLAT:this._segList.addPt(i.p1),this._segList.addPt(r.p1);break;case Cn.CAP_SQUARE:var u=new C;u.x=Math.abs(this._distance)*Math.cos(a),u.y=Math.abs(this._distance)*Math.sin(a);var l=new C(i.p1.x+u.x,i.p1.y+u.y),c=new C(r.p1.x+u.x,r.p1.y+u.y);this._segList.addPt(l),this._segList.addPt(c)}},Pn.prototype.getCoordinates=function(){return this._segList.getCoordinates()},Pn.prototype.addMitreJoin=function(t,e,n,i){var r=!0,o=null;try{o=k.intersection(e.p0,e.p1,n.p0,n.p1);(i<=0?1:o.distance(t)/Math.abs(i))>this._bufParams.getMitreLimit()&&(r=!1)}catch(t){if(!(t instanceof X))throw t;o=new C(0,0),r=!1}r?this._segList.addPt(o):this.addLimitedMitreJoin(e,n,i,this._bufParams.getMitreLimit())},Pn.prototype.addFilletCorner=function(t,e,n,i,r){var o=e.x-t.x,s=e.y-t.y,a=Math.atan2(s,o),u=n.x-t.x,l=n.y-t.y,c=Math.atan2(l,u);i===at.CLOCKWISE?a<=c&&(a+=2*Math.PI):a>=c&&(a-=2*Math.PI),this._segList.addPt(e),this.addFilletArc(t,a,c,i,r),this._segList.addPt(n)},Pn.prototype.addOutsideTurn=function(t,e){if(this._offset0.p1.distance(this._offset1.p0)<this._distance*Pn.OFFSET_SEGMENT_SEPARATION_FACTOR)return this._segList.addPt(this._offset0.p1),null;this._bufParams.getJoinStyle()===Cn.JOIN_MITRE?this.addMitreJoin(this._s1,this._offset0,this._offset1,this._distance):this._bufParams.getJoinStyle()===Cn.JOIN_BEVEL?this.addBevelJoin(this._offset0,this._offset1):(e&&this._segList.addPt(this._offset0.p1),this.addFilletCorner(this._s1,this._offset0.p1,this._offset1.p0,t,this._distance),this._segList.addPt(this._offset1.p0))},Pn.prototype.createSquare=function(t){this._segList.addPt(new C(t.x+this._distance,t.y+this._distance)),this._segList.addPt(new C(t.x+this._distance,t.y-this._distance)),this._segList.addPt(new C(t.x-this._distance,t.y-this._distance)),this._segList.addPt(new C(t.x-this._distance,t.y+this._distance)),this._segList.closeRing()},Pn.prototype.addSegments=function(t,e){this._segList.addPts(t,e)},Pn.prototype.addFirstSegment=function(){this._segList.addPt(this._offset1.p0)},Pn.prototype.addLastSegment=function(){this._segList.addPt(this._offset1.p1)},Pn.prototype.initSideSegments=function(t,e,n){this._s1=t,this._s2=e,this._side=n,this._seg1.setCoordinates(t,e),this.computeOffsetSegment(this._seg1,n,this._distance,this._offset1)},Pn.prototype.addLimitedMitreJoin=function(t,e,n,i){var r=this._seg0.p1,o=Tn.angle(r,this._seg0.p0),s=Tn.angleBetweenOriented(this._seg0.p0,r,this._seg1.p1)/2,a=Tn.normalize(o+s),u=Tn.normalize(a+Math.PI),l=i*n,c=n-l*Math.abs(Math.sin(s)),p=r.x+l*Math.cos(u),h=r.y+l*Math.sin(u),f=new C(p,h),g=new dn(r,f),d=g.pointAlongOffset(1,c),y=g.pointAlongOffset(1,-c);this._side===Se.LEFT?(this._segList.addPt(d),this._segList.addPt(y)):(this._segList.addPt(y),this._segList.addPt(d))},Pn.prototype.computeOffsetSegment=function(t,e,n,i){var r=e===Se.LEFT?1:-1,o=t.p1.x-t.p0.x,s=t.p1.y-t.p0.y,a=Math.sqrt(o*o+s*s),u=r*n*o/a,l=r*n*s/a;i.p0.x=t.p0.x-l,i.p0.y=t.p0.y+u,i.p1.x=t.p1.x-l,i.p1.y=t.p1.y+u},Pn.prototype.addFilletArc=function(t,e,n,i,r){var o=i===at.CLOCKWISE?-1:1,s=Math.abs(e-n),a=Math.trunc(s/this._filletAngleQuantum+.5);if(a<1)return null;for(var u=s/a,l=0,c=new C;l<s;){var p=e+o*l;c.x=t.x+r*Math.cos(p),c.y=t.y+r*Math.sin(p),this._segList.addPt(c),l+=u}},Pn.prototype.addInsideTurn=function(t,e){if(this._li.computeIntersection(this._offset0.p0,this._offset0.p1,this._offset1.p0,this._offset1.p1),this._li.hasIntersection())this._segList.addPt(this._li.getIntersection(0));else if(this._hasNarrowConcaveAngle=!0,this._offset0.p1.distance(this._offset1.p0)<this._distance*Pn.INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR)this._segList.addPt(this._offset0.p1);else{if(this._segList.addPt(this._offset0.p1),this._closingSegLengthFactor>0){var n=new C((this._closingSegLengthFactor*this._offset0.p1.x+this._s1.x)/(this._closingSegLengthFactor+1),(this._closingSegLengthFactor*this._offset0.p1.y+this._s1.y)/(this._closingSegLengthFactor+1));this._segList.addPt(n);var i=new C((this._closingSegLengthFactor*this._offset1.p0.x+this._s1.x)/(this._closingSegLengthFactor+1),(this._closingSegLengthFactor*this._offset1.p0.y+this._s1.y)/(this._closingSegLengthFactor+1));this._segList.addPt(i)}else this._segList.addPt(this._s1);this._segList.addPt(this._offset1.p0)}},Pn.prototype.createCircle=function(t){var e=new C(t.x+this._distance,t.y);this._segList.addPt(e),this.addFilletArc(t,0,2*Math.PI,-1,this._distance),this._segList.closeRing()},Pn.prototype.addBevelJoin=function(t,e){this._segList.addPt(t.p1),this._segList.addPt(e.p0)},Pn.prototype.init=function(t){this._distance=t,this._maxCurveSegmentError=t*(1-Math.cos(this._filletAngleQuantum/2)),this._segList=new wn,this._segList.setPrecisionModel(this._precisionModel),this._segList.setMinimumVertexDistance(t*Pn.CURVE_VERTEX_SNAP_DISTANCE_FACTOR)},Pn.prototype.addCollinear=function(t){this._li.computeIntersection(this._s0,this._s1,this._s1,this._s2);this._li.getIntersectionNum()>=2&&(this._bufParams.getJoinStyle()===Cn.JOIN_BEVEL||this._bufParams.getJoinStyle()===Cn.JOIN_MITRE?(t&&this._segList.addPt(this._offset0.p1),this._segList.addPt(this._offset1.p0)):this.addFilletCorner(this._s1,this._offset0.p1,this._offset1.p0,at.CLOCKWISE,this._distance))},Pn.prototype.closeRing=function(){this._segList.closeRing()},Pn.prototype.hasNarrowConcaveAngle=function(){return this._hasNarrowConcaveAngle},Pn.prototype.interfaces_=function(){return[]},Pn.prototype.getClass=function(){return Pn},Dn.OFFSET_SEGMENT_SEPARATION_FACTOR.get=function(){return.001},Dn.INSIDE_TURN_VERTEX_SNAP_DISTANCE_FACTOR.get=function(){return.001},Dn.CURVE_VERTEX_SNAP_DISTANCE_FACTOR.get=function(){return 1e-6},Dn.MAX_CLOSING_SEG_LEN_FACTOR.get=function(){return 80},Object.defineProperties(Pn,Dn);var Mn=function(){this._distance=0,this._precisionModel=null,this._bufParams=null;var t=arguments[0],e=arguments[1];this._precisionModel=t,this._bufParams=e};Mn.prototype.getOffsetCurve=function(t,e){if(this._distance=e,0===e)return null;var n=e<0,i=Math.abs(e),r=this.getSegGen(i);t.length<=1?this.computePointCurve(t[0],r):this.computeOffsetCurve(t,n,r);var o=r.getCoordinates();return n&&Lt.reverse(o),o},Mn.prototype.computeSingleSidedBufferCurve=function(t,e,n){var i=this.simplifyTolerance(this._distance);if(e){n.addSegments(t,!0);var r=Ln.simplify(t,-i),o=r.length-1;n.initSideSegments(r[o],r[o-1],Se.LEFT),n.addFirstSegment();for(var s=o-2;s>=0;s--)n.addNextSegment(r[s],!0)}else{n.addSegments(t,!1);var a=Ln.simplify(t,i),u=a.length-1;n.initSideSegments(a[0],a[1],Se.LEFT),n.addFirstSegment();for(var l=2;l<=u;l++)n.addNextSegment(a[l],!0)}n.addLastSegment(),n.closeRing()},Mn.prototype.computeRingBufferCurve=function(t,e,n){var i=this.simplifyTolerance(this._distance);e===Se.RIGHT&&(i=-i);var r=Ln.simplify(t,i),o=r.length-1;n.initSideSegments(r[o-1],r[0],e);for(var s=1;s<=o;s++){var a=1!==s;n.addNextSegment(r[s],a)}n.closeRing()},Mn.prototype.computeLineBufferCurve=function(t,e){var n=this.simplifyTolerance(this._distance),i=Ln.simplify(t,n),r=i.length-1;e.initSideSegments(i[0],i[1],Se.LEFT);for(var o=2;o<=r;o++)e.addNextSegment(i[o],!0);e.addLastSegment(),e.addLineEndCap(i[r-1],i[r]);var s=Ln.simplify(t,-n),a=s.length-1;e.initSideSegments(s[a],s[a-1],Se.LEFT);for(var u=a-2;u>=0;u--)e.addNextSegment(s[u],!0);e.addLastSegment(),e.addLineEndCap(s[1],s[0]),e.closeRing()},Mn.prototype.computePointCurve=function(t,e){switch(this._bufParams.getEndCapStyle()){case Cn.CAP_ROUND:e.createCircle(t);break;case Cn.CAP_SQUARE:e.createSquare(t)}},Mn.prototype.getLineCurve=function(t,e){if(this._distance=e,e<0&&!this._bufParams.isSingleSided())return null;if(0===e)return null;var n=Math.abs(e),i=this.getSegGen(n);if(t.length<=1)this.computePointCurve(t[0],i);else if(this._bufParams.isSingleSided()){var r=e<0;this.computeSingleSidedBufferCurve(t,r,i)}else this.computeLineBufferCurve(t,i);return i.getCoordinates()},Mn.prototype.getBufferParameters=function(){return this._bufParams},Mn.prototype.simplifyTolerance=function(t){return t*this._bufParams.getSimplifyFactor()},Mn.prototype.getRingCurve=function(t,e,n){if(this._distance=n,t.length<=2)return this.getLineCurve(t,n);if(0===n)return Mn.copyCoordinates(t);var i=this.getSegGen(n);return this.computeRingBufferCurve(t,e,i),i.getCoordinates()},Mn.prototype.computeOffsetCurve=function(t,e,n){var i=this.simplifyTolerance(this._distance);if(e){var r=Ln.simplify(t,-i),o=r.length-1;n.initSideSegments(r[o],r[o-1],Se.LEFT),n.addFirstSegment();for(var s=o-2;s>=0;s--)n.addNextSegment(r[s],!0)}else{var a=Ln.simplify(t,i),u=a.length-1;n.initSideSegments(a[0],a[1],Se.LEFT),n.addFirstSegment();for(var l=2;l<=u;l++)n.addNextSegment(a[l],!0)}n.addLastSegment()},Mn.prototype.getSegGen=function(t){return new Pn(this._precisionModel,this._bufParams,t)},Mn.prototype.interfaces_=function(){return[]},Mn.prototype.getClass=function(){return Mn},Mn.copyCoordinates=function(t){for(var e=new Array(t.length).fill(null),n=0;n<e.length;n++)e[n]=new C(t[n]);return e};var An=function(){this._subgraphs=null,this._seg=new dn,this._cga=new at;var t=arguments[0];this._subgraphs=t},Fn={DepthSegment:{configurable:!0}};An.prototype.findStabbedSegments=function(){if(1===arguments.length){for(var t=arguments[0],e=new Nt,n=this._subgraphs.iterator();n.hasNext();){var i=n.next(),r=i.getEnvelope();t.y<r.getMinY()||t.y>r.getMaxY()||this.findStabbedSegments(t,i.getDirectedEdges(),e)}return e}if(3===arguments.length)if(T(arguments[2],xt)&&arguments[0]instanceof C&&arguments[1]instanceof ze)for(var o=arguments[0],s=arguments[1],a=arguments[2],u=s.getEdge().getCoordinates(),l=0;l<u.length-1;l++){this._seg.p0=u[l],this._seg.p1=u[l+1],this._seg.p0.y>this._seg.p1.y&&this._seg.reverse();if(!(Math.max(this._seg.p0.x,this._seg.p1.x)<o.x)&&!(this._seg.isHorizontal()||o.y<this._seg.p0.y||o.y>this._seg.p1.y||at.computeOrientation(this._seg.p0,this._seg.p1,o)===at.RIGHT)){var c=s.getDepth(Se.LEFT);this._seg.p0.equals(u[l])||(c=s.getDepth(Se.RIGHT));var p=new Gn(this._seg,c);a.add(p)}}else if(T(arguments[2],xt)&&arguments[0]instanceof C&&T(arguments[1],xt))for(var h=arguments[0],f=arguments[1],g=arguments[2],d=f.iterator();d.hasNext();){var y=d.next();y.isForward()&&this.findStabbedSegments(h,y,g)}},An.prototype.getDepth=function(t){var e=this.findStabbedSegments(t);if(0===e.size())return 0;return $e.min(e)._leftDepth},An.prototype.interfaces_=function(){return[]},An.prototype.getClass=function(){return An},Fn.DepthSegment.get=function(){return Gn},Object.defineProperties(An,Fn);var Gn=function(){this._upwardSeg=null,this._leftDepth=null;var t=arguments[0],e=arguments[1];this._upwardSeg=new dn(t),this._leftDepth=e};Gn.prototype.compareTo=function(t){var e=t;if(this._upwardSeg.minX()>=e._upwardSeg.maxX())return 1;if(this._upwardSeg.maxX()<=e._upwardSeg.minX())return-1;var n=this._upwardSeg.orientationIndex(e._upwardSeg);return 0!==n?n:0!=(n=-1*e._upwardSeg.orientationIndex(this._upwardSeg))?n:this._upwardSeg.compareTo(e._upwardSeg)},Gn.prototype.compareX=function(t,e){var n=t.p0.compareTo(e.p0);return 0!==n?n:t.p1.compareTo(e.p1)},Gn.prototype.toString=function(){return this._upwardSeg.toString()},Gn.prototype.interfaces_=function(){return[E]},Gn.prototype.getClass=function(){return Gn};var qn=function(t,e,n){this.p0=t||null,this.p1=e||null,this.p2=n||null};qn.prototype.area=function(){return qn.area(this.p0,this.p1,this.p2)},qn.prototype.signedArea=function(){return qn.signedArea(this.p0,this.p1,this.p2)},qn.prototype.interpolateZ=function(t){if(null===t)throw new m("Supplied point is null.");return qn.interpolateZ(t,this.p0,this.p1,this.p2)},qn.prototype.longestSideLength=function(){return qn.longestSideLength(this.p0,this.p1,this.p2)},qn.prototype.isAcute=function(){return qn.isAcute(this.p0,this.p1,this.p2)},qn.prototype.circumcentre=function(){return qn.circumcentre(this.p0,this.p1,this.p2)},qn.prototype.area3D=function(){return qn.area3D(this.p0,this.p1,this.p2)},qn.prototype.centroid=function(){return qn.centroid(this.p0,this.p1,this.p2)},qn.prototype.inCentre=function(){return qn.inCentre(this.p0,this.p1,this.p2)},qn.prototype.interfaces_=function(){return[]},qn.prototype.getClass=function(){return qn},qn.area=function(t,e,n){return Math.abs(((n.x-t.x)*(e.y-t.y)-(e.x-t.x)*(n.y-t.y))/2)},qn.signedArea=function(t,e,n){return((n.x-t.x)*(e.y-t.y)-(e.x-t.x)*(n.y-t.y))/2},qn.det=function(t,e,n,i){return t*i-e*n},qn.interpolateZ=function(t,e,n,i){var r=e.x,o=e.y,s=n.x-r,a=i.x-r,u=n.y-o,l=i.y-o,c=s*l-a*u,p=t.x-r,h=t.y-o,f=(l*p-a*h)/c,g=(-u*p+s*h)/c;return e.z+f*(n.z-e.z)+g*(i.z-e.z)},qn.longestSideLength=function(t,e,n){var i=t.distance(e),r=e.distance(n),o=n.distance(t),s=i;return r>s&&(s=r),o>s&&(s=o),s},qn.isAcute=function(t,e,n){return!!Tn.isAcute(t,e,n)&&(!!Tn.isAcute(e,n,t)&&!!Tn.isAcute(n,t,e))},qn.circumcentre=function(t,e,n){var i=n.x,r=n.y,o=t.x-i,s=t.y-r,a=e.x-i,u=e.y-r,l=2*qn.det(o,s,a,u),c=qn.det(s,o*o+s*s,u,a*a+u*u),p=qn.det(o,o*o+s*s,a,a*a+u*u);return new C(i-c/l,r+p/l)},qn.perpendicularBisector=function(t,e){var n=e.x-t.x,i=e.y-t.y,r=new k(t.x+n/2,t.y+i/2,1),o=new k(t.x-i+n/2,t.y+n+i/2,1);return new k(r,o)},qn.angleBisector=function(t,e,n){var i=e.distance(t),r=i/(i+e.distance(n)),o=n.x-t.x,s=n.y-t.y;return new C(t.x+r*o,t.y+r*s)},qn.area3D=function(t,e,n){var i=e.x-t.x,r=e.y-t.y,o=e.z-t.z,s=n.x-t.x,a=n.y-t.y,u=n.z-t.z,l=r*u-o*a,c=o*s-i*u,p=i*a-r*s,h=l*l+c*c+p*p,f=Math.sqrt(h)/2;return f},qn.centroid=function(t,e,n){var i=(t.x+e.x+n.x)/3,r=(t.y+e.y+n.y)/3;return new C(i,r)},qn.inCentre=function(t,e,n){var i=e.distance(n),r=t.distance(n),o=t.distance(e),s=i+r+o,a=(i*t.x+r*e.x+o*n.x)/s,u=(i*t.y+r*e.y+o*n.y)/s;return new C(a,u)};var Bn=function(){this._inputGeom=null,this._distance=null,this._curveBuilder=null,this._curveList=new Nt;var t=arguments[0],e=arguments[1],n=arguments[2];this._inputGeom=t,this._distance=e,this._curveBuilder=n};Bn.prototype.addPoint=function(t){if(this._distance<=0)return null;var e=t.getCoordinates(),n=this._curveBuilder.getLineCurve(e,this._distance);this.addCurve(n,w.EXTERIOR,w.INTERIOR)},Bn.prototype.addPolygon=function(t){var e=this._distance,n=Se.LEFT;this._distance<0&&(e=-this._distance,n=Se.RIGHT);var i=t.getExteriorRing(),r=Lt.removeRepeatedPoints(i.getCoordinates());if(this._distance<0&&this.isErodedCompletely(i,this._distance))return null;if(this._distance<=0&&r.length<3)return null;this.addPolygonRing(r,e,n,w.EXTERIOR,w.INTERIOR);for(var o=0;o<t.getNumInteriorRing();o++){var s=t.getInteriorRingN(o),a=Lt.removeRepeatedPoints(s.getCoordinates());this._distance>0&&this.isErodedCompletely(s,-this._distance)||this.addPolygonRing(a,e,Se.opposite(n),w.INTERIOR,w.EXTERIOR)}},Bn.prototype.isTriangleErodedCompletely=function(t,e){var n=new qn(t[0],t[1],t[2]),i=n.inCentre();return at.distancePointLine(i,n.p0,n.p1)<Math.abs(e)},Bn.prototype.addLineString=function(t){if(this._distance<=0&&!this._curveBuilder.getBufferParameters().isSingleSided())return null;var e=Lt.removeRepeatedPoints(t.getCoordinates()),n=this._curveBuilder.getLineCurve(e,this._distance);this.addCurve(n,w.EXTERIOR,w.INTERIOR)},Bn.prototype.addCurve=function(t,e,n){if(null===t||t.length<2)return null;var i=new gn(t,new Pe(0,w.BOUNDARY,e,n));this._curveList.add(i)},Bn.prototype.getCurves=function(){return this.add(this._inputGeom),this._curveList},Bn.prototype.addPolygonRing=function(t,e,n,i,r){if(0===e&&t.length<ee.MINIMUM_VALID_SIZE)return null;var o=i,s=r;t.length>=ee.MINIMUM_VALID_SIZE&&at.isCCW(t)&&(o=r,s=i,n=Se.opposite(n));var a=this._curveBuilder.getRingCurve(t,n,e);this.addCurve(a,o,s)},Bn.prototype.add=function(t){if(t.isEmpty())return null;t instanceof $t?this.addPolygon(t):t instanceof Kt?this.addLineString(t):t instanceof Qt?this.addPoint(t):t instanceof te?this.addCollection(t):t instanceof Xt?this.addCollection(t):t instanceof ne?this.addCollection(t):t instanceof zt&&this.addCollection(t)},Bn.prototype.isErodedCompletely=function(t,e){var n=t.getCoordinates();if(n.length<4)return e<0;if(4===n.length)return this.isTriangleErodedCompletely(n,e);var i=t.getEnvelopeInternal(),r=Math.min(i.getHeight(),i.getWidth());return e<0&&2*Math.abs(e)>r},Bn.prototype.addCollection=function(t){for(var e=0;e<t.getNumGeometries();e++){var n=t.getGeometryN(e);this.add(n)}},Bn.prototype.interfaces_=function(){return[]},Bn.prototype.getClass=function(){return Bn};var Vn=function(){};Vn.prototype.locate=function(t){},Vn.prototype.interfaces_=function(){return[]},Vn.prototype.getClass=function(){return Vn};var Un=function(){this._parent=null,this._atStart=null,this._max=null,this._index=null,this._subcollectionIterator=null;var t=arguments[0];this._parent=t,this._atStart=!0,this._index=0,this._max=t.getNumGeometries()};Un.prototype.next=function(){if(this._atStart)return this._atStart=!1,Un.isAtomic(this._parent)&&this._index++,this._parent;if(null!==this._subcollectionIterator){if(this._subcollectionIterator.hasNext())return this._subcollectionIterator.next();this._subcollectionIterator=null}if(this._index>=this._max)throw new i;var t=this._parent.getGeometryN(this._index++);return t instanceof zt?(this._subcollectionIterator=new Un(t),this._subcollectionIterator.next()):t},Un.prototype.remove=function(){throw new Error(this.getClass().getName())},Un.prototype.hasNext=function(){if(this._atStart)return!0;if(null!==this._subcollectionIterator){if(this._subcollectionIterator.hasNext())return!0;this._subcollectionIterator=null}return!(this._index>=this._max)},Un.prototype.interfaces_=function(){return[Et]},Un.prototype.getClass=function(){return Un},Un.isAtomic=function(t){return!(t instanceof zt)};var zn=function(){this._geom=null;var t=arguments[0];this._geom=t};zn.prototype.locate=function(t){return zn.locate(t,this._geom)},zn.prototype.interfaces_=function(){return[Vn]},zn.prototype.getClass=function(){return zn},zn.isPointInRing=function(t,e){return!!e.getEnvelopeInternal().intersects(t)&&at.isPointInRing(t,e.getCoordinates())},zn.containsPointInPolygon=function(t,e){if(e.isEmpty())return!1;var n=e.getExteriorRing();if(!zn.isPointInRing(t,n))return!1;for(var i=0;i<e.getNumInteriorRing();i++){var r=e.getInteriorRingN(i);if(zn.isPointInRing(t,r))return!1}return!0},zn.containsPoint=function(t,e){if(e instanceof $t)return zn.containsPointInPolygon(t,e);if(e instanceof zt)for(var n=new Un(e);n.hasNext();){var i=n.next();if(i!==e&&zn.containsPoint(t,i))return!0}return!1},zn.locate=function(t,e){return e.isEmpty()?w.EXTERIOR:zn.containsPoint(t,e)?w.INTERIOR:w.EXTERIOR};var Xn=function(){this._edgeMap=new p,this._edgeList=null,this._ptInAreaLocation=[w.NONE,w.NONE]};Xn.prototype.getNextCW=function(t){this.getEdges();var e=this._edgeList.indexOf(t),n=e-1;return 0===e&&(n=this._edgeList.size()-1),this._edgeList.get(n)},Xn.prototype.propagateSideLabels=function(t){for(var e=w.NONE,n=this.iterator();n.hasNext();){var i=n.next().getLabel();i.isArea(t)&&i.getLocation(t,Se.LEFT)!==w.NONE&&(e=i.getLocation(t,Se.LEFT))}if(e===w.NONE)return null;for(var r=e,o=this.iterator();o.hasNext();){var s=o.next(),a=s.getLabel();if(a.getLocation(t,Se.ON)===w.NONE&&a.setLocation(t,Se.ON,r),a.isArea(t)){var u=a.getLocation(t,Se.LEFT),l=a.getLocation(t,Se.RIGHT);if(l!==w.NONE){if(l!==r)throw new we("side location conflict",s.getCoordinate());u===w.NONE&&et.shouldNeverReachHere("found single null side (at "+s.getCoordinate()+")"),r=u}else et.isTrue(a.getLocation(t,Se.LEFT)===w.NONE,"found single null side"),a.setLocation(t,Se.RIGHT,r),a.setLocation(t,Se.LEFT,r)}}},Xn.prototype.getCoordinate=function(){var t=this.iterator();if(!t.hasNext())return null;return t.next().getCoordinate()},Xn.prototype.print=function(t){Y.out.println("EdgeEndStar:   "+this.getCoordinate());for(var e=this.iterator();e.hasNext();){e.next().print(t)}},Xn.prototype.isAreaLabelsConsistent=function(t){return this.computeEdgeEndLabels(t.getBoundaryNodeRule()),this.checkAreaLabelsConsistent(0)},Xn.prototype.checkAreaLabelsConsistent=function(t){var e=this.getEdges();if(e.size()<=0)return!0;var n=e.size()-1,i=e.get(n).getLabel().getLocation(t,Se.LEFT);et.isTrue(i!==w.NONE,"Found unlabelled area edge");for(var r=i,o=this.iterator();o.hasNext();){var s=o.next().getLabel();et.isTrue(s.isArea(t),"Found non-area edge");var a=s.getLocation(t,Se.LEFT),u=s.getLocation(t,Se.RIGHT);if(a===u)return!1;if(u!==r)return!1;r=a}return!0},Xn.prototype.findIndex=function(t){this.iterator();for(var e=0;e<this._edgeList.size();e++){if(this._edgeList.get(e)===t)return e}return-1},Xn.prototype.iterator=function(){return this.getEdges().iterator()},Xn.prototype.getEdges=function(){return null===this._edgeList&&(this._edgeList=new Nt(this._edgeMap.values())),this._edgeList},Xn.prototype.getLocation=function(t,e,n){return this._ptInAreaLocation[t]===w.NONE&&(this._ptInAreaLocation[t]=zn.locate(e,n[t].getGeometry())),this._ptInAreaLocation[t]},Xn.prototype.toString=function(){var t=new D;t.append("EdgeEndStar:   "+this.getCoordinate()),t.append("\n");for(var e=this.iterator();e.hasNext();){var n=e.next();t.append(n),t.append("\n")}return t.toString()},Xn.prototype.computeEdgeEndLabels=function(t){for(var e=this.iterator();e.hasNext();){e.next().computeLabel(t)}},Xn.prototype.computeLabelling=function(t){this.computeEdgeEndLabels(t[0].getBoundaryNodeRule()),this.propagateSideLabels(0),this.propagateSideLabels(1);for(var e=[!1,!1],n=this.iterator();n.hasNext();)for(var i=n.next().getLabel(),r=0;r<2;r++)i.isLine(r)&&i.getLocation(r)===w.BOUNDARY&&(e[r]=!0);for(var o=this.iterator();o.hasNext();)for(var s=o.next(),a=s.getLabel(),u=0;u<2;u++)if(a.isAnyNull(u)){var l=w.NONE;if(e[u])l=w.EXTERIOR;else{var c=s.getCoordinate();l=this.getLocation(u,c,t)}a.setAllLocationsIfNull(u,l)}},Xn.prototype.getDegree=function(){return this._edgeMap.size()},Xn.prototype.insertEdgeEnd=function(t,e){this._edgeMap.put(t,e),this._edgeList=null},Xn.prototype.interfaces_=function(){return[]},Xn.prototype.getClass=function(){return Xn};var Yn=function(t){function e(){t.call(this),this._resultAreaEdgeList=null,this._label=null,this._SCANNING_FOR_INCOMING=1,this._LINKING_TO_OUTGOING=2}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.linkResultDirectedEdges=function(){this.getResultAreaEdges();for(var t=null,e=null,n=this._SCANNING_FOR_INCOMING,i=0;i<this._resultAreaEdgeList.size();i++){var r=this._resultAreaEdgeList.get(i),o=r.getSym();if(r.getLabel().isArea())switch(null===t&&r.isInResult()&&(t=r),n){case this._SCANNING_FOR_INCOMING:if(!o.isInResult())continue;e=o,n=this._LINKING_TO_OUTGOING;break;case this._LINKING_TO_OUTGOING:if(!r.isInResult())continue;e.setNext(r),n=this._SCANNING_FOR_INCOMING}}if(n===this._LINKING_TO_OUTGOING){if(null===t)throw new we("no outgoing dirEdge found",this.getCoordinate());et.isTrue(t.isInResult(),"unable to link last incoming dirEdge"),e.setNext(t)}},e.prototype.insert=function(t){var e=t;this.insertEdgeEnd(e,e)},e.prototype.getRightmostEdge=function(){var t=this.getEdges(),e=t.size();if(e<1)return null;var n=t.get(0);if(1===e)return n;var i=t.get(e-1),r=n.getQuadrant(),o=i.getQuadrant();return Be.isNorthern(r)&&Be.isNorthern(o)?n:Be.isNorthern(r)||Be.isNorthern(o)?0!==n.getDy()?n:0!==i.getDy()?i:(et.shouldNeverReachHere("found two horizontal edges incident on node"),null):i},e.prototype.print=function(t){Y.out.println("DirectedEdgeStar: "+this.getCoordinate());for(var e=this.iterator();e.hasNext();){var n=e.next();t.print("out "),n.print(t),t.println(),t.print("in "),n.getSym().print(t),t.println()}},e.prototype.getResultAreaEdges=function(){if(null!==this._resultAreaEdgeList)return this._resultAreaEdgeList;this._resultAreaEdgeList=new Nt;for(var t=this.iterator();t.hasNext();){var e=t.next();(e.isInResult()||e.getSym().isInResult())&&this._resultAreaEdgeList.add(e)}return this._resultAreaEdgeList},e.prototype.updateLabelling=function(t){for(var e=this.iterator();e.hasNext();){var n=e.next().getLabel();n.setAllLocationsIfNull(0,t.getLocation(0)),n.setAllLocationsIfNull(1,t.getLocation(1))}},e.prototype.linkAllDirectedEdges=function(){this.getEdges();for(var t=null,e=null,n=this._edgeList.size()-1;n>=0;n--){var i=this._edgeList.get(n),r=i.getSym();null===e&&(e=r),null!==t&&r.setNext(t),t=i}e.setNext(t)},e.prototype.computeDepths=function(){if(1===arguments.length){var t=arguments[0],e=this.findIndex(t),n=t.getDepth(Se.LEFT),i=t.getDepth(Se.RIGHT),r=this.computeDepths(e+1,this._edgeList.size(),n);if(this.computeDepths(0,e,r)!==i)throw new we("depth mismatch at "+t.getCoordinate())}else if(3===arguments.length){for(var o=arguments[0],s=arguments[1],a=arguments[2],u=o;u<s;u++){var l=this._edgeList.get(u);l.setEdgeDepths(Se.RIGHT,a),a=l.getDepth(Se.LEFT)}return a}},e.prototype.mergeSymLabels=function(){for(var t=this.iterator();t.hasNext();){var e=t.next();e.getLabel().merge(e.getSym().getLabel())}},e.prototype.linkMinimalDirectedEdges=function(t){for(var e=null,n=null,i=this._SCANNING_FOR_INCOMING,r=this._resultAreaEdgeList.size()-1;r>=0;r--){var o=this._resultAreaEdgeList.get(r),s=o.getSym();switch(null===e&&o.getEdgeRing()===t&&(e=o),i){case this._SCANNING_FOR_INCOMING:if(s.getEdgeRing()!==t)continue;n=s,i=this._LINKING_TO_OUTGOING;break;case this._LINKING_TO_OUTGOING:if(o.getEdgeRing()!==t)continue;n.setNextMin(o),i=this._SCANNING_FOR_INCOMING}}i===this._LINKING_TO_OUTGOING&&(et.isTrue(null!==e,"found null for first outgoing dirEdge"),et.isTrue(e.getEdgeRing()===t,"unable to link last incoming dirEdge"),n.setNextMin(e))},e.prototype.getOutgoingDegree=function(){if(0===arguments.length){for(var t=0,e=this.iterator();e.hasNext();){e.next().isInResult()&&t++}return t}if(1===arguments.length){for(var n=arguments[0],i=0,r=this.iterator();r.hasNext();){r.next().getEdgeRing()===n&&i++}return i}},e.prototype.getLabel=function(){return this._label},e.prototype.findCoveredLineEdges=function(){for(var t=w.NONE,e=this.iterator();e.hasNext();){var n=e.next(),i=n.getSym();if(!n.isLineEdge()){if(n.isInResult()){t=w.INTERIOR;break}if(i.isInResult()){t=w.EXTERIOR;break}}}if(t===w.NONE)return null;for(var r=t,o=this.iterator();o.hasNext();){var s=o.next(),a=s.getSym();s.isLineEdge()?s.getEdge().setCovered(r===w.INTERIOR):(s.isInResult()&&(r=w.EXTERIOR),a.isInResult()&&(r=w.INTERIOR))}},e.prototype.computeLabelling=function(e){t.prototype.computeLabelling.call(this,e),this._label=new Pe(w.NONE);for(var n=this.iterator();n.hasNext();)for(var i=n.next().getEdge().getLabel(),r=0;r<2;r++){var o=i.getLocation(r);o!==w.INTERIOR&&o!==w.BOUNDARY||this._label.setLocation(r,w.INTERIOR)}},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Xn),kn=function(t){function e(){t.apply(this,arguments)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.createNode=function(t){return new Ge(t,new Yn)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Xe),jn=function t(){this._pts=null,this._orientation=null;var e=arguments[0];this._pts=e,this._orientation=t.orientation(e)};jn.prototype.compareTo=function(t){var e=t;return jn.compareOriented(this._pts,this._orientation,e._pts,e._orientation)},jn.prototype.interfaces_=function(){return[E]},jn.prototype.getClass=function(){return jn},jn.orientation=function(t){return 1===Lt.increasingDirection(t)},jn.compareOriented=function(t,e,n,i){for(var r=e?1:-1,o=i?1:-1,s=e?t.length:-1,a=i?n.length:-1,u=e?0:t.length-1,l=i?0:n.length-1;;){var c=t[u].compareTo(n[l]);if(0!==c)return c;var p=(u+=r)===s,h=(l+=o)===a;if(p&&!h)return-1;if(!p&&h)return 1;if(p&&h)return 0}};var Hn=function(){this._edges=new Nt,this._ocaMap=new p};Hn.prototype.print=function(t){t.print("MULTILINESTRING ( ");for(var e=0;e<this._edges.size();e++){var n=this._edges.get(e);e>0&&t.print(","),t.print("(");for(var i=n.getCoordinates(),r=0;r<i.length;r++)r>0&&t.print(","),t.print(i[r].x+" "+i[r].y);t.println(")")}t.print(")  ")},Hn.prototype.addAll=function(t){for(var e=t.iterator();e.hasNext();)this.add(e.next())},Hn.prototype.findEdgeIndex=function(t){for(var e=0;e<this._edges.size();e++)if(this._edges.get(e).equals(t))return e;return-1},Hn.prototype.iterator=function(){return this._edges.iterator()},Hn.prototype.getEdges=function(){return this._edges},Hn.prototype.get=function(t){return this._edges.get(t)},Hn.prototype.findEqualEdge=function(t){var e=new jn(t.getCoordinates());return this._ocaMap.get(e)},Hn.prototype.add=function(t){this._edges.add(t);var e=new jn(t.getCoordinates());this._ocaMap.put(e,t)},Hn.prototype.interfaces_=function(){return[]},Hn.prototype.getClass=function(){return Hn};var Wn=function(){};Wn.prototype.processIntersections=function(t,e,n,i){},Wn.prototype.isDone=function(){},Wn.prototype.interfaces_=function(){return[]},Wn.prototype.getClass=function(){return Wn};var Kn=function(){this._hasIntersection=!1,this._hasProper=!1,this._hasProperInterior=!1,this._hasInterior=!1,this._properIntersectionPoint=null,this._li=null,this._isSelfIntersection=null,this.numIntersections=0,this.numInteriorIntersections=0,this.numProperIntersections=0,this.numTests=0;var t=arguments[0];this._li=t};Kn.prototype.isTrivialIntersection=function(t,e,n,i){if(t===n&&1===this._li.getIntersectionNum()){if(Kn.isAdjacentSegments(e,i))return!0;if(t.isClosed()){var r=t.size()-1;if(0===e&&i===r||0===i&&e===r)return!0}}return!1},Kn.prototype.getProperIntersectionPoint=function(){return this._properIntersectionPoint},Kn.prototype.hasProperInteriorIntersection=function(){return this._hasProperInterior},Kn.prototype.getLineIntersector=function(){return this._li},Kn.prototype.hasProperIntersection=function(){return this._hasProper},Kn.prototype.processIntersections=function(t,e,n,i){if(t===n&&e===i)return null;this.numTests++;var r=t.getCoordinates()[e],o=t.getCoordinates()[e+1],s=n.getCoordinates()[i],a=n.getCoordinates()[i+1];this._li.computeIntersection(r,o,s,a),this._li.hasIntersection()&&(this.numIntersections++,this._li.isInteriorIntersection()&&(this.numInteriorIntersections++,this._hasInterior=!0),this.isTrivialIntersection(t,e,n,i)||(this._hasIntersection=!0,t.addIntersections(this._li,e,0),n.addIntersections(this._li,i,1),this._li.isProper()&&(this.numProperIntersections++,this._hasProper=!0,this._hasProperInterior=!0)))},Kn.prototype.hasIntersection=function(){return this._hasIntersection},Kn.prototype.isDone=function(){return!1},Kn.prototype.hasInteriorIntersection=function(){return this._hasInterior},Kn.prototype.interfaces_=function(){return[Wn]},Kn.prototype.getClass=function(){return Kn},Kn.isAdjacentSegments=function(t,e){return 1===Math.abs(t-e)};var Jn=function(){this.coord=null,this.segmentIndex=null,this.dist=null;var t=arguments[0],e=arguments[1],n=arguments[2];this.coord=new C(t),this.segmentIndex=e,this.dist=n};Jn.prototype.getSegmentIndex=function(){return this.segmentIndex},Jn.prototype.getCoordinate=function(){return this.coord},Jn.prototype.print=function(t){t.print(this.coord),t.print(" seg # = "+this.segmentIndex),t.println(" dist = "+this.dist)},Jn.prototype.compareTo=function(t){var e=t;return this.compare(e.segmentIndex,e.dist)},Jn.prototype.isEndPoint=function(t){return 0===this.segmentIndex&&0===this.dist||this.segmentIndex===t},Jn.prototype.toString=function(){return this.coord+" seg # = "+this.segmentIndex+" dist = "+this.dist},Jn.prototype.getDistance=function(){return this.dist},Jn.prototype.compare=function(t,e){return this.segmentIndex<t?-1:this.segmentIndex>t?1:this.dist<e?-1:this.dist>e?1:0},Jn.prototype.interfaces_=function(){return[E]},Jn.prototype.getClass=function(){return Jn};var Qn=function(){this._nodeMap=new p,this.edge=null;var t=arguments[0];this.edge=t};Qn.prototype.print=function(t){t.println("Intersections:");for(var e=this.iterator();e.hasNext();){e.next().print(t)}},Qn.prototype.iterator=function(){return this._nodeMap.values().iterator()},Qn.prototype.addSplitEdges=function(t){this.addEndpoints();for(var e=this.iterator(),n=e.next();e.hasNext();){var i=e.next(),r=this.createSplitEdge(n,i);t.add(r),n=i}},Qn.prototype.addEndpoints=function(){var t=this.edge.pts.length-1;this.add(this.edge.pts[0],0,0),this.add(this.edge.pts[t],t,0)},Qn.prototype.createSplitEdge=function(t,e){var n=e.segmentIndex-t.segmentIndex+2,i=this.edge.pts[e.segmentIndex],r=e.dist>0||!e.coord.equals2D(i);r||n--;var o=new Array(n).fill(null),s=0;o[s++]=new C(t.coord);for(var a=t.segmentIndex+1;a<=e.segmentIndex;a++)o[s++]=this.edge.pts[a];return r&&(o[s]=e.coord),new ni(o,new Pe(this.edge._label))},Qn.prototype.add=function(t,e,n){var i=new Jn(t,e,n),r=this._nodeMap.get(i);return null!==r?r:(this._nodeMap.put(i,i),i)},Qn.prototype.isIntersection=function(t){for(var e=this.iterator();e.hasNext();){if(e.next().coord.equals(t))return!0}return!1},Qn.prototype.interfaces_=function(){return[]},Qn.prototype.getClass=function(){return Qn};var Zn=function(){};Zn.prototype.getChainStartIndices=function(t){var e=0,n=new Nt;n.add(new M(e));do{var i=this.findChainEnd(t,e);n.add(new M(i)),e=i}while(e<t.length-1);return Zn.toIntArray(n)},Zn.prototype.findChainEnd=function(t,e){for(var n=Be.quadrant(t[e],t[e+1]),i=e+1;i<t.length;){if(Be.quadrant(t[i-1],t[i])!==n)break;i++}return i-1},Zn.prototype.interfaces_=function(){return[]},Zn.prototype.getClass=function(){return Zn},Zn.toIntArray=function(t){for(var e=new Array(t.size()).fill(null),n=0;n<e.length;n++)e[n]=t.get(n).intValue();return e};var $n=function(){this.e=null,this.pts=null,this.startIndex=null,this.env1=new j,this.env2=new j;var t=arguments[0];this.e=t,this.pts=t.getCoordinates();var e=new Zn;this.startIndex=e.getChainStartIndices(this.pts)};$n.prototype.getCoordinates=function(){return this.pts},$n.prototype.getMaxX=function(t){var e=this.pts[this.startIndex[t]].x,n=this.pts[this.startIndex[t+1]].x;return e>n?e:n},$n.prototype.getMinX=function(t){var e=this.pts[this.startIndex[t]].x,n=this.pts[this.startIndex[t+1]].x;return e<n?e:n},$n.prototype.computeIntersectsForChain=function(){if(4===arguments.length){var t=arguments[0],e=arguments[1],n=arguments[2],i=arguments[3];this.computeIntersectsForChain(this.startIndex[t],this.startIndex[t+1],e,e.startIndex[n],e.startIndex[n+1],i)}else if(6===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2],a=arguments[3],u=arguments[4],l=arguments[5],c=this.pts[r],p=this.pts[o],h=s.pts[a],f=s.pts[u];if(o-r==1&&u-a==1)return l.addIntersections(this.e,r,s.e,a),null;if(this.env1.init(c,p),this.env2.init(h,f),!this.env1.intersects(this.env2))return null;var g=Math.trunc((r+o)/2),d=Math.trunc((a+u)/2);r<g&&(a<d&&this.computeIntersectsForChain(r,g,s,a,d,l),d<u&&this.computeIntersectsForChain(r,g,s,d,u,l)),g<o&&(a<d&&this.computeIntersectsForChain(g,o,s,a,d,l),d<u&&this.computeIntersectsForChain(g,o,s,d,u,l))}},$n.prototype.getStartIndexes=function(){return this.startIndex},$n.prototype.computeIntersects=function(t,e){for(var n=0;n<this.startIndex.length-1;n++)for(var i=0;i<t.startIndex.length-1;i++)this.computeIntersectsForChain(n,t,i,e)},$n.prototype.interfaces_=function(){return[]},$n.prototype.getClass=function(){return $n};var ti=function t(){this._depth=Array(2).fill().map(function(){return Array(3)});for(var e=0;e<2;e++)for(var n=0;n<3;n++)this._depth[e][n]=t.NULL_VALUE},ei={NULL_VALUE:{configurable:!0}};ti.prototype.getDepth=function(t,e){return this._depth[t][e]},ti.prototype.setDepth=function(t,e,n){this._depth[t][e]=n},ti.prototype.isNull=function(){if(0===arguments.length){for(var t=0;t<2;t++)for(var e=0;e<3;e++)if(this._depth[t][e]!==ti.NULL_VALUE)return!1;return!0}if(1===arguments.length){var n=arguments[0];return this._depth[n][1]===ti.NULL_VALUE}if(2===arguments.length){var i=arguments[0],r=arguments[1];return this._depth[i][r]===ti.NULL_VALUE}},ti.prototype.normalize=function(){for(var t=0;t<2;t++)if(!this.isNull(t)){var e=this._depth[t][1];this._depth[t][2]<e&&(e=this._depth[t][2]),e<0&&(e=0);for(var n=1;n<3;n++){var i=0;this._depth[t][n]>e&&(i=1),this._depth[t][n]=i}}},ti.prototype.getDelta=function(t){return this._depth[t][Se.RIGHT]-this._depth[t][Se.LEFT]},ti.prototype.getLocation=function(t,e){return this._depth[t][e]<=0?w.EXTERIOR:w.INTERIOR},ti.prototype.toString=function(){return"A: "+this._depth[0][1]+","+this._depth[0][2]+" B: "+this._depth[1][1]+","+this._depth[1][2]},ti.prototype.add=function(){if(1===arguments.length)for(var t=arguments[0],e=0;e<2;e++)for(var n=1;n<3;n++){var i=t.getLocation(e,n);i!==w.EXTERIOR&&i!==w.INTERIOR||(this.isNull(e,n)?this._depth[e][n]=ti.depthAtLocation(i):this._depth[e][n]+=ti.depthAtLocation(i))}else if(3===arguments.length){var r=arguments[0],o=arguments[1];arguments[2]===w.INTERIOR&&this._depth[r][o]++}},ti.prototype.interfaces_=function(){return[]},ti.prototype.getClass=function(){return ti},ti.depthAtLocation=function(t){return t===w.EXTERIOR?0:t===w.INTERIOR?1:ti.NULL_VALUE},ei.NULL_VALUE.get=function(){return-1},Object.defineProperties(ti,ei);var ni=function(t){function e(){if(t.call(this),this.pts=null,this._env=null,this.eiList=new Qn(this),this._name=null,this._mce=null,this._isIsolated=!0,this._depth=new ti,this._depthDelta=0,1===arguments.length){var n=arguments[0];e.call(this,n,null)}else if(2===arguments.length){var i=arguments[0],r=arguments[1];this.pts=i,this._label=r}}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.getDepth=function(){return this._depth},e.prototype.getCollapsedEdge=function(){var t=new Array(2).fill(null);t[0]=this.pts[0],t[1]=this.pts[1];return new e(t,Pe.toLineLabel(this._label))},e.prototype.isIsolated=function(){return this._isIsolated},e.prototype.getCoordinates=function(){return this.pts},e.prototype.setIsolated=function(t){this._isIsolated=t},e.prototype.setName=function(t){this._name=t},e.prototype.equals=function(t){if(!(t instanceof e))return!1;var n=t;if(this.pts.length!==n.pts.length)return!1;for(var i=!0,r=!0,o=this.pts.length,s=0;s<this.pts.length;s++)if(this.pts[s].equals2D(n.pts[s])||(i=!1),this.pts[s].equals2D(n.pts[--o])||(r=!1),!i&&!r)return!1;return!0},e.prototype.getCoordinate=function(){if(0===arguments.length)return this.pts.length>0?this.pts[0]:null;if(1===arguments.length){var t=arguments[0];return this.pts[t]}},e.prototype.print=function(t){t.print("edge "+this._name+": "),t.print("LINESTRING (");for(var e=0;e<this.pts.length;e++)e>0&&t.print(","),t.print(this.pts[e].x+" "+this.pts[e].y);t.print(")  "+this._label+" "+this._depthDelta)},e.prototype.computeIM=function(t){e.updateIM(this._label,t)},e.prototype.isCollapsed=function(){return!!this._label.isArea()&&(3===this.pts.length&&!!this.pts[0].equals(this.pts[2]))},e.prototype.isClosed=function(){return this.pts[0].equals(this.pts[this.pts.length-1])},e.prototype.getMaximumSegmentIndex=function(){return this.pts.length-1},e.prototype.getDepthDelta=function(){return this._depthDelta},e.prototype.getNumPoints=function(){return this.pts.length},e.prototype.printReverse=function(t){t.print("edge "+this._name+": ");for(var e=this.pts.length-1;e>=0;e--)t.print(this.pts[e]+" ");t.println("")},e.prototype.getMonotoneChainEdge=function(){return null===this._mce&&(this._mce=new $n(this)),this._mce},e.prototype.getEnvelope=function(){if(null===this._env){this._env=new j;for(var t=0;t<this.pts.length;t++)this._env.expandToInclude(this.pts[t])}return this._env},e.prototype.addIntersection=function(t,e,n,i){var r=new C(t.getIntersection(i)),o=e,s=t.getEdgeDistance(n,i),a=o+1;if(a<this.pts.length){var u=this.pts[a];r.equals2D(u)&&(o=a,s=0)}this.eiList.add(r,o,s)},e.prototype.toString=function(){var t=new D;t.append("edge "+this._name+": "),t.append("LINESTRING (");for(var e=0;e<this.pts.length;e++)e>0&&t.append(","),t.append(this.pts[e].x+" "+this.pts[e].y);return t.append(")  "+this._label+" "+this._depthDelta),t.toString()},e.prototype.isPointwiseEqual=function(t){if(this.pts.length!==t.pts.length)return!1;for(var e=0;e<this.pts.length;e++)if(!this.pts[e].equals2D(t.pts[e]))return!1;return!0},e.prototype.setDepthDelta=function(t){this._depthDelta=t},e.prototype.getEdgeIntersectionList=function(){return this.eiList},e.prototype.addIntersections=function(t,e,n){for(var i=0;i<t.getIntersectionNum();i++)this.addIntersection(t,e,n,i)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e.updateIM=function(){if(2!==arguments.length)return t.prototype.updateIM.apply(this,arguments);var e=arguments[0],n=arguments[1];n.setAtLeastIfValid(e.getLocation(0,Se.ON),e.getLocation(1,Se.ON),1),e.isArea()&&(n.setAtLeastIfValid(e.getLocation(0,Se.LEFT),e.getLocation(1,Se.LEFT),2),n.setAtLeastIfValid(e.getLocation(0,Se.RIGHT),e.getLocation(1,Se.RIGHT),2))},e}(Fe),ii=function(t){this._workingPrecisionModel=null,this._workingNoder=null,this._geomFact=null,this._graph=null,this._edgeList=new Hn,this._bufParams=t||null};ii.prototype.setWorkingPrecisionModel=function(t){this._workingPrecisionModel=t},ii.prototype.insertUniqueEdge=function(t){var e=this._edgeList.findEqualEdge(t);if(null!==e){var n=e.getLabel(),i=t.getLabel();e.isPointwiseEqual(t)||(i=new Pe(t.getLabel())).flip(),n.merge(i);var r=ii.depthDelta(i),o=e.getDepthDelta()+r;e.setDepthDelta(o)}else this._edgeList.add(t),t.setDepthDelta(ii.depthDelta(t.getLabel()))},ii.prototype.buildSubgraphs=function(t,e){for(var n=new Nt,i=t.iterator();i.hasNext();){var r=i.next(),o=r.getRightmostCoordinate(),s=new An(n).getDepth(o);r.computeDepth(s),r.findResultEdges(),n.add(r),e.add(r.getDirectedEdges(),r.getNodes())}},ii.prototype.createSubgraphs=function(t){for(var e=new Nt,n=t.getNodes().iterator();n.hasNext();){var i=n.next();if(!i.isVisited()){var r=new Te;r.create(i),e.add(r)}}return $e.sort(e,$e.reverseOrder()),e},ii.prototype.createEmptyResultGeometry=function(){return this._geomFact.createPolygon()},ii.prototype.getNoder=function(t){if(null!==this._workingNoder)return this._workingNoder;var e=new xn,n=new rt;return n.setPrecisionModel(t),e.setSegmentIntersector(new Kn(n)),e},ii.prototype.buffer=function(t,e){var n=this._workingPrecisionModel;null===n&&(n=t.getPrecisionModel()),this._geomFact=t.getFactory();var i=new Mn(n,this._bufParams),r=new Bn(t,e,i).getCurves();if(r.size()<=0)return this.createEmptyResultGeometry();this.computeNodedEdges(r,n),this._graph=new Ye(new kn),this._graph.addEdges(this._edgeList.getEdges());var o=this.createSubgraphs(this._graph),s=new ke(this._geomFact);this.buildSubgraphs(o,s);var a=s.getPolygons();if(a.size()<=0)return this.createEmptyResultGeometry();return this._geomFact.buildGeometry(a)},ii.prototype.computeNodedEdges=function(t,e){var n=this.getNoder(e);n.computeNodes(t);for(var i=n.getNodedSubstrings().iterator();i.hasNext();){var r=i.next(),o=r.getCoordinates();if(2!==o.length||!o[0].equals2D(o[1])){var s=r.getData(),a=new ni(r.getCoordinates(),new Pe(s));this.insertUniqueEdge(a)}}},ii.prototype.setNoder=function(t){this._workingNoder=t},ii.prototype.interfaces_=function(){return[]},ii.prototype.getClass=function(){return ii},ii.depthDelta=function(t){var e=t.getLocation(0,Se.LEFT),n=t.getLocation(0,Se.RIGHT);return e===w.INTERIOR&&n===w.EXTERIOR?1:e===w.EXTERIOR&&n===w.INTERIOR?-1:0},ii.convertSegStrings=function(t){for(var e=new _e,n=new Nt;t.hasNext();){var i=t.next(),r=e.createLineString(i.getCoordinates());n.add(r)}return e.buildGeometry(n)};var ri=function(){if(this._noder=null,this._scaleFactor=null,this._offsetX=null,this._offsetY=null,this._isScaled=!1,2===arguments.length){var t=arguments[0],e=arguments[1];this._noder=t,this._scaleFactor=e,this._offsetX=0,this._offsetY=0,this._isScaled=!this.isIntegerPrecision()}else if(4===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2],o=arguments[3];this._noder=n,this._scaleFactor=i,this._offsetX=r,this._offsetY=o,this._isScaled=!this.isIntegerPrecision()}};ri.prototype.rescale=function(){if(T(arguments[0],It))for(var t=arguments[0].iterator();t.hasNext();){var e=t.next();this.rescale(e.getCoordinates())}else if(arguments[0]instanceof Array){for(var n=arguments[0],i=0;i<n.length;i++)n[i].x=n[i].x/this._scaleFactor+this._offsetX,n[i].y=n[i].y/this._scaleFactor+this._offsetY;2===n.length&&n[0].equals2D(n[1])&&Y.out.println(n)}},ri.prototype.scale=function(){if(T(arguments[0],It)){for(var t=arguments[0],e=new Nt,n=t.iterator();n.hasNext();){var i=n.next();e.add(new gn(this.scale(i.getCoordinates()),i.getData()))}return e}if(arguments[0]instanceof Array){for(var r=arguments[0],o=new Array(r.length).fill(null),s=0;s<r.length;s++)o[s]=new C(Math.round((r[s].x-this._offsetX)*this._scaleFactor),Math.round((r[s].y-this._offsetY)*this._scaleFactor),r[s].z);return Lt.removeRepeatedPoints(o)}},ri.prototype.isIntegerPrecision=function(){return 1===this._scaleFactor},ri.prototype.getNodedSubstrings=function(){var t=this._noder.getNodedSubstrings();return this._isScaled&&this.rescale(t),t},ri.prototype.computeNodes=function(t){var e=t;this._isScaled&&(e=this.scale(t)),this._noder.computeNodes(e)},ri.prototype.interfaces_=function(){return[In]},ri.prototype.getClass=function(){return ri};var oi=function(){this._li=new rt,this._segStrings=null;var t=arguments[0];this._segStrings=t},si={fact:{configurable:!0}};oi.prototype.checkEndPtVertexIntersections=function(){if(0===arguments.length)for(var t=this._segStrings.iterator();t.hasNext();){var e=t.next().getCoordinates();this.checkEndPtVertexIntersections(e[0],this._segStrings),this.checkEndPtVertexIntersections(e[e.length-1],this._segStrings)}else if(2===arguments.length)for(var n=arguments[0],i=arguments[1].iterator();i.hasNext();)for(var r=i.next().getCoordinates(),o=1;o<r.length-1;o++)if(r[o].equals(n))throw new $("found endpt/interior pt intersection at index "+o+" :pt "+n)},oi.prototype.checkInteriorIntersections=function(){if(0===arguments.length)for(var t=this._segStrings.iterator();t.hasNext();)for(var e=t.next(),n=this._segStrings.iterator();n.hasNext();){var i=n.next();this.checkInteriorIntersections(e,i)}else if(2===arguments.length)for(var r=arguments[0],o=arguments[1],s=r.getCoordinates(),a=o.getCoordinates(),u=0;u<s.length-1;u++)for(var l=0;l<a.length-1;l++)this.checkInteriorIntersections(r,u,o,l);else if(4===arguments.length){var c=arguments[0],p=arguments[1],h=arguments[2],f=arguments[3];if(c===h&&p===f)return null;var g=c.getCoordinates()[p],d=c.getCoordinates()[p+1],y=h.getCoordinates()[f],_=h.getCoordinates()[f+1];if(this._li.computeIntersection(g,d,y,_),this._li.hasIntersection()&&(this._li.isProper()||this.hasInteriorIntersection(this._li,g,d)||this.hasInteriorIntersection(this._li,y,_)))throw new $("found non-noded intersection at "+g+"-"+d+" and "+y+"-"+_)}},oi.prototype.checkValid=function(){this.checkEndPtVertexIntersections(),this.checkInteriorIntersections(),this.checkCollapses()},oi.prototype.checkCollapses=function(){if(0===arguments.length)for(var t=this._segStrings.iterator();t.hasNext();){var e=t.next();this.checkCollapses(e)}else if(1===arguments.length)for(var n=arguments[0].getCoordinates(),i=0;i<n.length-2;i++)this.checkCollapse(n[i],n[i+1],n[i+2])},oi.prototype.hasInteriorIntersection=function(t,e,n){for(var i=0;i<t.getIntersectionNum();i++){var r=t.getIntersection(i);if(!r.equals(e)&&!r.equals(n))return!0}return!1},oi.prototype.checkCollapse=function(t,e,n){if(t.equals(n))throw new $("found non-noded collapse at "+oi.fact.createLineString([t,e,n]))},oi.prototype.interfaces_=function(){return[]},oi.prototype.getClass=function(){return oi},si.fact.get=function(){return new _e},Object.defineProperties(oi,si);var ai=function(){this._li=null,this._pt=null,this._originalPt=null,this._ptScaled=null,this._p0Scaled=null,this._p1Scaled=null,this._scaleFactor=null,this._minx=null,this._maxx=null,this._miny=null,this._maxy=null,this._corner=new Array(4).fill(null),this._safeEnv=null;var t=arguments[0],e=arguments[1],n=arguments[2];if(this._originalPt=t,this._pt=t,this._scaleFactor=e,this._li=n,e<=0)throw new m("Scale factor must be non-zero");1!==e&&(this._pt=new C(this.scale(t.x),this.scale(t.y)),this._p0Scaled=new C,this._p1Scaled=new C),this.initCorners(this._pt)},ui={SAFE_ENV_EXPANSION_FACTOR:{configurable:!0}};ai.prototype.intersectsScaled=function(t,e){var n=Math.min(t.x,e.x),i=Math.max(t.x,e.x),r=Math.min(t.y,e.y),o=Math.max(t.y,e.y),s=this._maxx<n||this._minx>i||this._maxy<r||this._miny>o;if(s)return!1;var a=this.intersectsToleranceSquare(t,e);return et.isTrue(!(s&&a),"Found bad envelope test"),a},ai.prototype.initCorners=function(t){this._minx=t.x-.5,this._maxx=t.x+.5,this._miny=t.y-.5,this._maxy=t.y+.5,this._corner[0]=new C(this._maxx,this._maxy),this._corner[1]=new C(this._minx,this._maxy),this._corner[2]=new C(this._minx,this._miny),this._corner[3]=new C(this._maxx,this._miny)},ai.prototype.intersects=function(t,e){return 1===this._scaleFactor?this.intersectsScaled(t,e):(this.copyScaled(t,this._p0Scaled),this.copyScaled(e,this._p1Scaled),this.intersectsScaled(this._p0Scaled,this._p1Scaled))},ai.prototype.scale=function(t){return Math.round(t*this._scaleFactor)},ai.prototype.getCoordinate=function(){return this._originalPt},ai.prototype.copyScaled=function(t,e){e.x=this.scale(t.x),e.y=this.scale(t.y)},ai.prototype.getSafeEnvelope=function(){if(null===this._safeEnv){var t=ai.SAFE_ENV_EXPANSION_FACTOR/this._scaleFactor;this._safeEnv=new j(this._originalPt.x-t,this._originalPt.x+t,this._originalPt.y-t,this._originalPt.y+t)}return this._safeEnv},ai.prototype.intersectsPixelClosure=function(t,e){return this._li.computeIntersection(t,e,this._corner[0],this._corner[1]),!!this._li.hasIntersection()||(this._li.computeIntersection(t,e,this._corner[1],this._corner[2]),!!this._li.hasIntersection()||(this._li.computeIntersection(t,e,this._corner[2],this._corner[3]),!!this._li.hasIntersection()||(this._li.computeIntersection(t,e,this._corner[3],this._corner[0]),!!this._li.hasIntersection())))},ai.prototype.intersectsToleranceSquare=function(t,e){var n=!1,i=!1;return this._li.computeIntersection(t,e,this._corner[0],this._corner[1]),!!this._li.isProper()||(this._li.computeIntersection(t,e,this._corner[1],this._corner[2]),!!this._li.isProper()||(this._li.hasIntersection()&&(n=!0),this._li.computeIntersection(t,e,this._corner[2],this._corner[3]),!!this._li.isProper()||(this._li.hasIntersection()&&(i=!0),this._li.computeIntersection(t,e,this._corner[3],this._corner[0]),!!this._li.isProper()||(!(!n||!i)||(!!t.equals(this._pt)||!!e.equals(this._pt))))))},ai.prototype.addSnappedNode=function(t,e){var n=t.getCoordinate(e),i=t.getCoordinate(e+1);return!!this.intersects(n,i)&&(t.addIntersection(this.getCoordinate(),e),!0)},ai.prototype.interfaces_=function(){return[]},ai.prototype.getClass=function(){return ai},ui.SAFE_ENV_EXPANSION_FACTOR.get=function(){return.75},Object.defineProperties(ai,ui);var li=function(){this.tempEnv1=new j,this.selectedSegment=new dn};li.prototype.select=function(){if(1===arguments.length);else if(2===arguments.length){var t=arguments[0],e=arguments[1];t.getLineSegment(e,this.selectedSegment),this.select(this.selectedSegment)}},li.prototype.interfaces_=function(){return[]},li.prototype.getClass=function(){return li};var ci=function(){this._index=null;var t=arguments[0];this._index=t},pi={HotPixelSnapAction:{configurable:!0}};ci.prototype.snap=function(){if(1===arguments.length){var t=arguments[0];return this.snap(t,null,-1)}if(3===arguments.length){var e=arguments[0],n=arguments[1],i=arguments[2],r=e.getSafeEnvelope(),o=new hi(e,n,i);return this._index.query(r,{interfaces_:function(){return[Ke]},visitItem:function(t){t.select(r,o)}}),o.isNodeAdded()}},ci.prototype.interfaces_=function(){return[]},ci.prototype.getClass=function(){return ci},pi.HotPixelSnapAction.get=function(){return hi},Object.defineProperties(ci,pi);var hi=function(t){function e(){t.call(this),this._hotPixel=null,this._parentEdge=null,this._hotPixelVertexIndex=null,this._isNodeAdded=!1;var e=arguments[0],n=arguments[1],i=arguments[2];this._hotPixel=e,this._parentEdge=n,this._hotPixelVertexIndex=i}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.isNodeAdded=function(){return this._isNodeAdded},e.prototype.select=function(){if(2!==arguments.length)return t.prototype.select.apply(this,arguments);var e=arguments[0],n=arguments[1],i=e.getContext();if(null!==this._parentEdge&&i===this._parentEdge&&n===this._hotPixelVertexIndex)return null;this._isNodeAdded=this._hotPixel.addSnappedNode(i,n)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(li),fi=function(){this._li=null,this._interiorIntersections=null;var t=arguments[0];this._li=t,this._interiorIntersections=new Nt};fi.prototype.processIntersections=function(t,e,n,i){if(t===n&&e===i)return null;var r=t.getCoordinates()[e],o=t.getCoordinates()[e+1],s=n.getCoordinates()[i],a=n.getCoordinates()[i+1];if(this._li.computeIntersection(r,o,s,a),this._li.hasIntersection()&&this._li.isInteriorIntersection()){for(var u=0;u<this._li.getIntersectionNum();u++)this._interiorIntersections.add(this._li.getIntersection(u));t.addIntersections(this._li,e,0),n.addIntersections(this._li,i,1)}},fi.prototype.isDone=function(){return!1},fi.prototype.getInteriorIntersections=function(){return this._interiorIntersections},fi.prototype.interfaces_=function(){return[Wn]},fi.prototype.getClass=function(){return fi};var gi=function(){this._pm=null,this._li=null,this._scaleFactor=null,this._noder=null,this._pointSnapper=null,this._nodedSegStrings=null;var t=arguments[0];this._pm=t,this._li=new rt,this._li.setPrecisionModel(t),this._scaleFactor=t.getScale()};gi.prototype.checkCorrectness=function(t){var e=gn.getNodedSubstrings(t),n=new oi(e);try{n.checkValid()}catch(t){if(!(t instanceof z))throw t;t.printStackTrace()}},gi.prototype.getNodedSubstrings=function(){return gn.getNodedSubstrings(this._nodedSegStrings)},gi.prototype.snapRound=function(t,e){var n=this.findInteriorIntersections(t,e);this.computeIntersectionSnaps(n),this.computeVertexSnaps(t)},gi.prototype.findInteriorIntersections=function(t,e){var n=new fi(e);return this._noder.setSegmentIntersector(n),this._noder.computeNodes(t),n.getInteriorIntersections()},gi.prototype.computeVertexSnaps=function(){if(T(arguments[0],It))for(var t=arguments[0].iterator();t.hasNext();){var e=t.next();this.computeVertexSnaps(e)}else if(arguments[0]instanceof gn)for(var n=arguments[0],i=n.getCoordinates(),r=0;r<i.length;r++){var o=new ai(i[r],this._scaleFactor,this._li);this._pointSnapper.snap(o,n,r)&&n.addIntersection(i[r],r)}},gi.prototype.computeNodes=function(t){this._nodedSegStrings=t,this._noder=new xn,this._pointSnapper=new ci(this._noder.getIndex()),this.snapRound(t,this._li)},gi.prototype.computeIntersectionSnaps=function(t){for(var e=t.iterator();e.hasNext();){var n=e.next(),i=new ai(n,this._scaleFactor,this._li);this._pointSnapper.snap(i)}},gi.prototype.interfaces_=function(){return[In]},gi.prototype.getClass=function(){return gi};var di=function(){if(this._argGeom=null,this._distance=null,this._bufParams=new Cn,this._resultGeometry=null,this._saveException=null,1===arguments.length){var t=arguments[0];this._argGeom=t}else if(2===arguments.length){var e=arguments[0],n=arguments[1];this._argGeom=e,this._bufParams=n}},yi={CAP_ROUND:{configurable:!0},CAP_BUTT:{configurable:!0},CAP_FLAT:{configurable:!0},CAP_SQUARE:{configurable:!0},MAX_PRECISION_DIGITS:{configurable:!0}};di.prototype.bufferFixedPrecision=function(t){var e=new ri(new gi(new fe(1)),t.getScale()),n=new ii(this._bufParams);n.setWorkingPrecisionModel(t),n.setNoder(e),this._resultGeometry=n.buffer(this._argGeom,this._distance)},di.prototype.bufferReducedPrecision=function(){var t=this;if(0===arguments.length){for(var e=di.MAX_PRECISION_DIGITS;e>=0;e--){try{t.bufferReducedPrecision(e)}catch(e){if(!(e instanceof we))throw e;t._saveException=e}if(null!==t._resultGeometry)return null}throw this._saveException}if(1===arguments.length){var n=arguments[0],i=di.precisionScaleFactor(this._argGeom,this._distance,n),r=new fe(i);this.bufferFixedPrecision(r)}},di.prototype.computeGeometry=function(){if(this.bufferOriginalPrecision(),null!==this._resultGeometry)return null;var t=this._argGeom.getFactory().getPrecisionModel();t.getType()===fe.FIXED?this.bufferFixedPrecision(t):this.bufferReducedPrecision()},di.prototype.setQuadrantSegments=function(t){this._bufParams.setQuadrantSegments(t)},di.prototype.bufferOriginalPrecision=function(){try{var t=new ii(this._bufParams);this._resultGeometry=t.buffer(this._argGeom,this._distance)}catch(t){if(!(t instanceof $))throw t;this._saveException=t}},di.prototype.getResultGeometry=function(t){return this._distance=t,this.computeGeometry(),this._resultGeometry},di.prototype.setEndCapStyle=function(t){this._bufParams.setEndCapStyle(t)},di.prototype.interfaces_=function(){return[]},di.prototype.getClass=function(){return di},di.bufferOp=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];return new di(t).getResultGeometry(e)}if(3===arguments.length){if(Number.isInteger(arguments[2])&&arguments[0]instanceof ct&&"number"==typeof arguments[1]){var n=arguments[0],i=arguments[1],r=arguments[2],o=new di(n);o.setQuadrantSegments(r);return o.getResultGeometry(i)}if(arguments[2]instanceof Cn&&arguments[0]instanceof ct&&"number"==typeof arguments[1]){var s=arguments[0],a=arguments[1],u=arguments[2];return new di(s,u).getResultGeometry(a)}}else if(4===arguments.length){var l=arguments[0],c=arguments[1],p=arguments[2],h=arguments[3],f=new di(l);f.setQuadrantSegments(p),f.setEndCapStyle(h);return f.getResultGeometry(c)}},di.precisionScaleFactor=function(t,e,n){var i=t.getEnvelopeInternal(),r=R.max(Math.abs(i.getMaxX()),Math.abs(i.getMaxY()),Math.abs(i.getMinX()),Math.abs(i.getMinY()))+2*(e>0?e:0),o=n-Math.trunc(Math.log(r)/Math.log(10)+1);return Math.pow(10,o)},yi.CAP_ROUND.get=function(){return Cn.CAP_ROUND},yi.CAP_BUTT.get=function(){return Cn.CAP_FLAT},yi.CAP_FLAT.get=function(){return Cn.CAP_FLAT},yi.CAP_SQUARE.get=function(){return Cn.CAP_SQUARE},yi.MAX_PRECISION_DIGITS.get=function(){return 12},Object.defineProperties(di,yi);var _i=function(){this._pt=[new C,new C],this._distance=v.NaN,this._isNull=!0};_i.prototype.getCoordinates=function(){return this._pt},_i.prototype.getCoordinate=function(t){return this._pt[t]},_i.prototype.setMinimum=function(){if(1===arguments.length){var t=arguments[0];this.setMinimum(t._pt[0],t._pt[1])}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(this._isNull)return this.initialize(e,n),null;var i=e.distance(n);i<this._distance&&this.initialize(e,n,i)}},_i.prototype.initialize=function(){if(0===arguments.length)this._isNull=!0;else if(2===arguments.length){var t=arguments[0],e=arguments[1];this._pt[0].setCoordinate(t),this._pt[1].setCoordinate(e),this._distance=t.distance(e),this._isNull=!1}else if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2];this._pt[0].setCoordinate(n),this._pt[1].setCoordinate(i),this._distance=r,this._isNull=!1}},_i.prototype.getDistance=function(){return this._distance},_i.prototype.setMaximum=function(){if(1===arguments.length){var t=arguments[0];this.setMaximum(t._pt[0],t._pt[1])}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(this._isNull)return this.initialize(e,n),null;var i=e.distance(n);i>this._distance&&this.initialize(e,n,i)}},_i.prototype.interfaces_=function(){return[]},_i.prototype.getClass=function(){return _i};var mi=function(){};mi.prototype.interfaces_=function(){return[]},mi.prototype.getClass=function(){return mi},mi.computeDistance=function(){if(arguments[2]instanceof _i&&arguments[0]instanceof Kt&&arguments[1]instanceof C)for(var t=arguments[0],e=arguments[1],n=arguments[2],i=t.getCoordinates(),r=new dn,o=0;o<i.length-1;o++){r.setCoordinates(i[o],i[o+1]);var s=r.closestPoint(e);n.setMinimum(s,e)}else if(arguments[2]instanceof _i&&arguments[0]instanceof $t&&arguments[1]instanceof C){var a=arguments[0],u=arguments[1],l=arguments[2];mi.computeDistance(a.getExteriorRing(),u,l);for(var c=0;c<a.getNumInteriorRing();c++)mi.computeDistance(a.getInteriorRingN(c),u,l)}else if(arguments[2]instanceof _i&&arguments[0]instanceof ct&&arguments[1]instanceof C){var p=arguments[0],h=arguments[1],f=arguments[2];if(p instanceof Kt)mi.computeDistance(p,h,f);else if(p instanceof $t)mi.computeDistance(p,h,f);else if(p instanceof zt)for(var g=p,d=0;d<g.getNumGeometries();d++){var y=g.getGeometryN(d);mi.computeDistance(y,h,f)}else f.setMinimum(p.getCoordinate(),h)}else if(arguments[2]instanceof _i&&arguments[0]instanceof dn&&arguments[1]instanceof C){var _=arguments[0],m=arguments[1],v=arguments[2],I=_.closestPoint(m);v.setMinimum(I,m)}};var vi=function(t){this._maxPtDist=new _i,this._inputGeom=t||null},Ii={MaxPointDistanceFilter:{configurable:!0},MaxMidpointDistanceFilter:{configurable:!0}};vi.prototype.computeMaxMidpointDistance=function(t){var e=new xi(this._inputGeom);t.apply(e),this._maxPtDist.setMaximum(e.getMaxPointDistance())},vi.prototype.computeMaxVertexDistance=function(t){var e=new Ei(this._inputGeom);t.apply(e),this._maxPtDist.setMaximum(e.getMaxPointDistance())},vi.prototype.findDistance=function(t){return this.computeMaxVertexDistance(t),this.computeMaxMidpointDistance(t),this._maxPtDist.getDistance()},vi.prototype.getDistancePoints=function(){return this._maxPtDist},vi.prototype.interfaces_=function(){return[]},vi.prototype.getClass=function(){return vi},Ii.MaxPointDistanceFilter.get=function(){return Ei},Ii.MaxMidpointDistanceFilter.get=function(){return xi},Object.defineProperties(vi,Ii);var Ei=function(t){this._maxPtDist=new _i,this._minPtDist=new _i,this._geom=t||null};Ei.prototype.filter=function(t){this._minPtDist.initialize(),mi.computeDistance(this._geom,t,this._minPtDist),this._maxPtDist.setMaximum(this._minPtDist)},Ei.prototype.getMaxPointDistance=function(){return this._maxPtDist},Ei.prototype.interfaces_=function(){return[ft]},Ei.prototype.getClass=function(){return Ei};var xi=function(t){this._maxPtDist=new _i,this._minPtDist=new _i,this._geom=t||null};xi.prototype.filter=function(t,e){if(0===e)return null;var n=t.getCoordinate(e-1),i=t.getCoordinate(e),r=new C((n.x+i.x)/2,(n.y+i.y)/2);this._minPtDist.initialize(),mi.computeDistance(this._geom,r,this._minPtDist),this._maxPtDist.setMaximum(this._minPtDist)},xi.prototype.isDone=function(){return!1},xi.prototype.isGeometryChanged=function(){return!1},xi.prototype.getMaxPointDistance=function(){return this._maxPtDist},xi.prototype.interfaces_=function(){return[Ut]},xi.prototype.getClass=function(){return xi};var Ni=function(t){this._comps=t||null};Ni.prototype.filter=function(t){t instanceof $t&&this._comps.add(t)},Ni.prototype.interfaces_=function(){return[Vt]},Ni.prototype.getClass=function(){return Ni},Ni.getPolygons=function(){if(1===arguments.length){var t=arguments[0];return Ni.getPolygons(t,new Nt)}if(2===arguments.length){var e=arguments[0],n=arguments[1];return e instanceof $t?n.add(e):e instanceof zt&&e.apply(new Ni(n)),n}};var Ci=function(){if(this._lines=null,this._isForcedToLineString=!1,1===arguments.length){var t=arguments[0];this._lines=t}else if(2===arguments.length){var e=arguments[0],n=arguments[1];this._lines=e,this._isForcedToLineString=n}};Ci.prototype.filter=function(t){if(this._isForcedToLineString&&t instanceof ee){var e=t.getFactory().createLineString(t.getCoordinateSequence());return this._lines.add(e),null}t instanceof Kt&&this._lines.add(t)},Ci.prototype.setForceToLineString=function(t){this._isForcedToLineString=t},Ci.prototype.interfaces_=function(){return[lt]},Ci.prototype.getClass=function(){return Ci},Ci.getGeometry=function(){if(1===arguments.length){var t=arguments[0];return t.getFactory().buildGeometry(Ci.getLines(t))}if(2===arguments.length){var e=arguments[0],n=arguments[1];return e.getFactory().buildGeometry(Ci.getLines(e,n))}},Ci.getLines=function(){if(1===arguments.length){var t=arguments[0];return Ci.getLines(t,!1)}if(2===arguments.length){if(T(arguments[0],It)&&T(arguments[1],It)){for(var e=arguments[0],n=arguments[1],i=e.iterator();i.hasNext();){var r=i.next();Ci.getLines(r,n)}return n}if(arguments[0]instanceof ct&&"boolean"==typeof arguments[1]){var o=arguments[0],s=arguments[1],a=new Nt;return o.apply(new Ci(a,s)),a}if(arguments[0]instanceof ct&&T(arguments[1],It)){var u=arguments[0],l=arguments[1];return u instanceof Kt?l.add(u):u.apply(new Ci(l)),l}}else if(3===arguments.length){if("boolean"==typeof arguments[2]&&T(arguments[0],It)&&T(arguments[1],It)){for(var c=arguments[0],p=arguments[1],h=arguments[2],f=c.iterator();f.hasNext();){var g=f.next();Ci.getLines(g,p,h)}return p}if("boolean"==typeof arguments[2]&&arguments[0]instanceof ct&&T(arguments[1],It)){var d=arguments[0],y=arguments[1],_=arguments[2];return d.apply(new Ci(y,_)),y}}};var Si=function(){if(this._boundaryRule=gt.OGC_SFS_BOUNDARY_RULE,this._isIn=null,this._numBoundaries=null,0===arguments.length);else if(1===arguments.length){var t=arguments[0];if(null===t)throw new m("Rule must be non-null");this._boundaryRule=t}};Si.prototype.locateInternal=function(){if(arguments[0]instanceof C&&arguments[1]instanceof $t){var t=arguments[0],e=arguments[1];if(e.isEmpty())return w.EXTERIOR;var n=e.getExteriorRing(),i=this.locateInPolygonRing(t,n);if(i===w.EXTERIOR)return w.EXTERIOR;if(i===w.BOUNDARY)return w.BOUNDARY;for(var r=0;r<e.getNumInteriorRing();r++){var o=e.getInteriorRingN(r),s=this.locateInPolygonRing(t,o);if(s===w.INTERIOR)return w.EXTERIOR;if(s===w.BOUNDARY)return w.BOUNDARY}return w.INTERIOR}if(arguments[0]instanceof C&&arguments[1]instanceof Kt){var a=arguments[0],u=arguments[1];if(!u.getEnvelopeInternal().intersects(a))return w.EXTERIOR;var l=u.getCoordinates();return u.isClosed()||!a.equals(l[0])&&!a.equals(l[l.length-1])?at.isOnLine(a,l)?w.INTERIOR:w.EXTERIOR:w.BOUNDARY}if(arguments[0]instanceof C&&arguments[1]instanceof Qt){var c=arguments[0];return arguments[1].getCoordinate().equals2D(c)?w.INTERIOR:w.EXTERIOR}},Si.prototype.locateInPolygonRing=function(t,e){return e.getEnvelopeInternal().intersects(t)?at.locatePointInRing(t,e.getCoordinates()):w.EXTERIOR},Si.prototype.intersects=function(t,e){return this.locate(t,e)!==w.EXTERIOR},Si.prototype.updateLocationInfo=function(t){t===w.INTERIOR&&(this._isIn=!0),t===w.BOUNDARY&&this._numBoundaries++},Si.prototype.computeLocation=function(t,e){if(e instanceof Qt&&this.updateLocationInfo(this.locateInternal(t,e)),e instanceof Kt)this.updateLocationInfo(this.locateInternal(t,e));else if(e instanceof $t)this.updateLocationInfo(this.locateInternal(t,e));else if(e instanceof Xt)for(var n=e,i=0;i<n.getNumGeometries();i++){var r=n.getGeometryN(i);this.updateLocationInfo(this.locateInternal(t,r))}else if(e instanceof ne)for(var o=e,s=0;s<o.getNumGeometries();s++){var a=o.getGeometryN(s);this.updateLocationInfo(this.locateInternal(t,a))}else if(e instanceof zt)for(var u=new Un(e);u.hasNext();){var l=u.next();l!==e&&this.computeLocation(t,l)}},Si.prototype.locate=function(t,e){return e.isEmpty()?w.EXTERIOR:e instanceof Kt?this.locateInternal(t,e):e instanceof $t?this.locateInternal(t,e):(this._isIn=!1,this._numBoundaries=0,this.computeLocation(t,e),this._boundaryRule.isInBoundary(this._numBoundaries)?w.BOUNDARY:this._numBoundaries>0||this._isIn?w.INTERIOR:w.EXTERIOR)},Si.prototype.interfaces_=function(){return[]},Si.prototype.getClass=function(){return Si};var Li=function t(){if(this._component=null,this._segIndex=null,this._pt=null,2===arguments.length){var e=arguments[0],n=arguments[1];t.call(this,e,t.INSIDE_AREA,n)}else if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2];this._component=i,this._segIndex=r,this._pt=o}},bi={INSIDE_AREA:{configurable:!0}};Li.prototype.isInsideArea=function(){return this._segIndex===Li.INSIDE_AREA},Li.prototype.getCoordinate=function(){return this._pt},Li.prototype.getGeometryComponent=function(){return this._component},Li.prototype.getSegmentIndex=function(){return this._segIndex},Li.prototype.interfaces_=function(){return[]},Li.prototype.getClass=function(){return Li},bi.INSIDE_AREA.get=function(){return-1},Object.defineProperties(Li,bi);var wi=function(t){this._pts=t||null};wi.prototype.filter=function(t){t instanceof Qt&&this._pts.add(t)},wi.prototype.interfaces_=function(){return[Vt]},wi.prototype.getClass=function(){return wi},wi.getPoints=function(){if(1===arguments.length){var t=arguments[0];return t instanceof Qt?$e.singletonList(t):wi.getPoints(t,new Nt)}if(2===arguments.length){var e=arguments[0],n=arguments[1];return e instanceof Qt?n.add(e):e instanceof zt&&e.apply(new wi(n)),n}};var Oi=function(){this._locations=null;var t=arguments[0];this._locations=t};Oi.prototype.filter=function(t){(t instanceof Qt||t instanceof Kt||t instanceof $t)&&this._locations.add(new Li(t,0,t.getCoordinate()))},Oi.prototype.interfaces_=function(){return[Vt]},Oi.prototype.getClass=function(){return Oi},Oi.getLocations=function(t){var e=new Nt;return t.apply(new Oi(e)),e};var Ti=function(){if(this._geom=null,this._terminateDistance=0,this._ptLocator=new Si,this._minDistanceLocation=null,this._minDistance=v.MAX_VALUE,2===arguments.length){var t=arguments[0],e=arguments[1];this._geom=[t,e],this._terminateDistance=0}else if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2];this._geom=new Array(2).fill(null),this._geom[0]=n,this._geom[1]=i,this._terminateDistance=r}};Ti.prototype.computeContainmentDistance=function(){if(0===arguments.length){var t=new Array(2).fill(null);if(this.computeContainmentDistance(0,t),this._minDistance<=this._terminateDistance)return null;this.computeContainmentDistance(1,t)}else if(2===arguments.length){var e=arguments[0],n=arguments[1],i=1-e,r=Ni.getPolygons(this._geom[e]);if(r.size()>0){var o=Oi.getLocations(this._geom[i]);if(this.computeContainmentDistance(o,r,n),this._minDistance<=this._terminateDistance)return this._minDistanceLocation[i]=n[0],this._minDistanceLocation[e]=n[1],null}}else if(3===arguments.length)if(arguments[2]instanceof Array&&T(arguments[0],xt)&&T(arguments[1],xt)){for(var s=arguments[0],a=arguments[1],u=arguments[2],l=0;l<s.size();l++)for(var c=s.get(l),p=0;p<a.size();p++)if(this.computeContainmentDistance(c,a.get(p),u),this._minDistance<=this._terminateDistance)return null}else if(arguments[2]instanceof Array&&arguments[0]instanceof Li&&arguments[1]instanceof $t){var h=arguments[0],f=arguments[1],g=arguments[2],d=h.getCoordinate();if(w.EXTERIOR!==this._ptLocator.locate(d,f))return this._minDistance=0,g[0]=h,g[1]=new Li(f,d),null}},Ti.prototype.computeMinDistanceLinesPoints=function(t,e,n){for(var i=0;i<t.size();i++)for(var r=t.get(i),o=0;o<e.size();o++){var s=e.get(o);if(this.computeMinDistance(r,s,n),this._minDistance<=this._terminateDistance)return null}},Ti.prototype.computeFacetDistance=function(){var t=new Array(2).fill(null),e=Ci.getLines(this._geom[0]),n=Ci.getLines(this._geom[1]),i=wi.getPoints(this._geom[0]),r=wi.getPoints(this._geom[1]);return this.computeMinDistanceLines(e,n,t),this.updateMinDistance(t,!1),this._minDistance<=this._terminateDistance?null:(t[0]=null,t[1]=null,this.computeMinDistanceLinesPoints(e,r,t),this.updateMinDistance(t,!1),this._minDistance<=this._terminateDistance?null:(t[0]=null,t[1]=null,this.computeMinDistanceLinesPoints(n,i,t),this.updateMinDistance(t,!0),this._minDistance<=this._terminateDistance?null:(t[0]=null,t[1]=null,this.computeMinDistancePoints(i,r,t),void this.updateMinDistance(t,!1))))},Ti.prototype.nearestLocations=function(){return this.computeMinDistance(),this._minDistanceLocation},Ti.prototype.updateMinDistance=function(t,e){if(null===t[0])return null;e?(this._minDistanceLocation[0]=t[1],this._minDistanceLocation[1]=t[0]):(this._minDistanceLocation[0]=t[0],this._minDistanceLocation[1]=t[1])},Ti.prototype.nearestPoints=function(){this.computeMinDistance();return[this._minDistanceLocation[0].getCoordinate(),this._minDistanceLocation[1].getCoordinate()]},Ti.prototype.computeMinDistance=function(){if(0===arguments.length){if(null!==this._minDistanceLocation)return null;if(this._minDistanceLocation=new Array(2).fill(null),this.computeContainmentDistance(),this._minDistance<=this._terminateDistance)return null;this.computeFacetDistance()}else if(3===arguments.length)if(arguments[2]instanceof Array&&arguments[0]instanceof Kt&&arguments[1]instanceof Qt){var t=arguments[0],e=arguments[1],n=arguments[2];if(t.getEnvelopeInternal().distance(e.getEnvelopeInternal())>this._minDistance)return null;for(var i=t.getCoordinates(),r=e.getCoordinate(),o=0;o<i.length-1;o++){var s=at.distancePointLine(r,i[o],i[o+1]);if(s<this._minDistance){this._minDistance=s;var a=new dn(i[o],i[o+1]).closestPoint(r);n[0]=new Li(t,o,a),n[1]=new Li(e,0,r)}if(this._minDistance<=this._terminateDistance)return null}}else if(arguments[2]instanceof Array&&arguments[0]instanceof Kt&&arguments[1]instanceof Kt){var u=arguments[0],l=arguments[1],c=arguments[2];if(u.getEnvelopeInternal().distance(l.getEnvelopeInternal())>this._minDistance)return null;for(var p=u.getCoordinates(),h=l.getCoordinates(),f=0;f<p.length-1;f++)for(var g=0;g<h.length-1;g++){var d=at.distanceLineLine(p[f],p[f+1],h[g],h[g+1]);if(d<this._minDistance){this._minDistance=d;var y=new dn(p[f],p[f+1]),_=new dn(h[g],h[g+1]),m=y.closestPoints(_);c[0]=new Li(u,f,m[0]),c[1]=new Li(l,g,m[1])}if(this._minDistance<=this._terminateDistance)return null}}},Ti.prototype.computeMinDistancePoints=function(t,e,n){for(var i=0;i<t.size();i++)for(var r=t.get(i),o=0;o<e.size();o++){var s=e.get(o),a=r.getCoordinate().distance(s.getCoordinate());if(a<this._minDistance&&(this._minDistance=a,n[0]=new Li(r,0,r.getCoordinate()),n[1]=new Li(s,0,s.getCoordinate())),this._minDistance<=this._terminateDistance)return null}},Ti.prototype.distance=function(){if(null===this._geom[0]||null===this._geom[1])throw new m("null geometries are not supported");return this._geom[0].isEmpty()||this._geom[1].isEmpty()?0:(this.computeMinDistance(),this._minDistance)},Ti.prototype.computeMinDistanceLines=function(t,e,n){for(var i=0;i<t.size();i++)for(var r=t.get(i),o=0;o<e.size();o++){var s=e.get(o);if(this.computeMinDistance(r,s,n),this._minDistance<=this._terminateDistance)return null}},Ti.prototype.interfaces_=function(){return[]},Ti.prototype.getClass=function(){return Ti},Ti.distance=function(t,e){return new Ti(t,e).distance()},Ti.isWithinDistance=function(t,e,n){return new Ti(t,e,n).distance()<=n},Ti.nearestPoints=function(t,e){return new Ti(t,e).nearestPoints()};var Ri=function(){this._pt=[new C,new C],this._distance=v.NaN,this._isNull=!0};Ri.prototype.getCoordinates=function(){return this._pt},Ri.prototype.getCoordinate=function(t){return this._pt[t]},Ri.prototype.setMinimum=function(){if(1===arguments.length){var t=arguments[0];this.setMinimum(t._pt[0],t._pt[1])}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(this._isNull)return this.initialize(e,n),null;var i=e.distance(n);i<this._distance&&this.initialize(e,n,i)}},Ri.prototype.initialize=function(){if(0===arguments.length)this._isNull=!0;else if(2===arguments.length){var t=arguments[0],e=arguments[1];this._pt[0].setCoordinate(t),this._pt[1].setCoordinate(e),this._distance=t.distance(e),this._isNull=!1}else if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2];this._pt[0].setCoordinate(n),this._pt[1].setCoordinate(i),this._distance=r,this._isNull=!1}},Ri.prototype.toString=function(){return Z.toLineString(this._pt[0],this._pt[1])},Ri.prototype.getDistance=function(){return this._distance},Ri.prototype.setMaximum=function(){if(1===arguments.length){var t=arguments[0];this.setMaximum(t._pt[0],t._pt[1])}else if(2===arguments.length){var e=arguments[0],n=arguments[1];if(this._isNull)return this.initialize(e,n),null;var i=e.distance(n);i>this._distance&&this.initialize(e,n,i)}},Ri.prototype.interfaces_=function(){return[]},Ri.prototype.getClass=function(){return Ri};var Pi=function(){};Pi.prototype.interfaces_=function(){return[]},Pi.prototype.getClass=function(){return Pi},Pi.computeDistance=function(){if(arguments[2]instanceof Ri&&arguments[0]instanceof Kt&&arguments[1]instanceof C)for(var t=arguments[0],e=arguments[1],n=arguments[2],i=new dn,r=t.getCoordinates(),o=0;o<r.length-1;o++){i.setCoordinates(r[o],r[o+1]);var s=i.closestPoint(e);n.setMinimum(s,e)}else if(arguments[2]instanceof Ri&&arguments[0]instanceof $t&&arguments[1]instanceof C){var a=arguments[0],u=arguments[1],l=arguments[2];Pi.computeDistance(a.getExteriorRing(),u,l);for(var c=0;c<a.getNumInteriorRing();c++)Pi.computeDistance(a.getInteriorRingN(c),u,l)}else if(arguments[2]instanceof Ri&&arguments[0]instanceof ct&&arguments[1]instanceof C){var p=arguments[0],h=arguments[1],f=arguments[2];if(p instanceof Kt)Pi.computeDistance(p,h,f);else if(p instanceof $t)Pi.computeDistance(p,h,f);else if(p instanceof zt)for(var g=p,d=0;d<g.getNumGeometries();d++){var y=g.getGeometryN(d);Pi.computeDistance(y,h,f)}else f.setMinimum(p.getCoordinate(),h)}else if(arguments[2]instanceof Ri&&arguments[0]instanceof dn&&arguments[1]instanceof C){var _=arguments[0],m=arguments[1],v=arguments[2],I=_.closestPoint(m);v.setMinimum(I,m)}};var Di=function(){this._g0=null,this._g1=null,this._ptDist=new Ri,this._densifyFrac=0;var t=arguments[0],e=arguments[1];this._g0=t,this._g1=e},Mi={MaxPointDistanceFilter:{configurable:!0},MaxDensifiedByFractionDistanceFilter:{configurable:!0}};Di.prototype.getCoordinates=function(){return this._ptDist.getCoordinates()},Di.prototype.setDensifyFraction=function(t){if(t>1||t<=0)throw new m("Fraction is not in range (0.0 - 1.0]");this._densifyFrac=t},Di.prototype.compute=function(t,e){this.computeOrientedDistance(t,e,this._ptDist),this.computeOrientedDistance(e,t,this._ptDist)},Di.prototype.distance=function(){return this.compute(this._g0,this._g1),this._ptDist.getDistance()},Di.prototype.computeOrientedDistance=function(t,e,n){var i=new Ai(e);if(t.apply(i),n.setMaximum(i.getMaxPointDistance()),this._densifyFrac>0){var r=new Fi(e,this._densifyFrac);t.apply(r),n.setMaximum(r.getMaxPointDistance())}},Di.prototype.orientedDistance=function(){return this.computeOrientedDistance(this._g0,this._g1,this._ptDist),this._ptDist.getDistance()},Di.prototype.interfaces_=function(){return[]},Di.prototype.getClass=function(){return Di},Di.distance=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];return new Di(t,e).distance()}if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2],o=new Di(n,i);return o.setDensifyFraction(r),o.distance()}},Mi.MaxPointDistanceFilter.get=function(){return Ai},Mi.MaxDensifiedByFractionDistanceFilter.get=function(){return Fi},Object.defineProperties(Di,Mi);var Ai=function(){this._maxPtDist=new Ri,this._minPtDist=new Ri,this._euclideanDist=new Pi,this._geom=null;var t=arguments[0];this._geom=t};Ai.prototype.filter=function(t){this._minPtDist.initialize(),Pi.computeDistance(this._geom,t,this._minPtDist),this._maxPtDist.setMaximum(this._minPtDist)},Ai.prototype.getMaxPointDistance=function(){return this._maxPtDist},Ai.prototype.interfaces_=function(){return[ft]},Ai.prototype.getClass=function(){return Ai};var Fi=function(){this._maxPtDist=new Ri,this._minPtDist=new Ri,this._geom=null,this._numSubSegs=0;var t=arguments[0],e=arguments[1];this._geom=t,this._numSubSegs=Math.trunc(Math.round(1/e))};Fi.prototype.filter=function(t,e){if(0===e)return null;for(var n=t.getCoordinate(e-1),i=t.getCoordinate(e),r=(i.x-n.x)/this._numSubSegs,o=(i.y-n.y)/this._numSubSegs,s=0;s<this._numSubSegs;s++){var a=n.x+s*r,u=n.y+s*o,l=new C(a,u);this._minPtDist.initialize(),Pi.computeDistance(this._geom,l,this._minPtDist),this._maxPtDist.setMaximum(this._minPtDist)}},Fi.prototype.isDone=function(){return!1},Fi.prototype.isGeometryChanged=function(){return!1},Fi.prototype.getMaxPointDistance=function(){return this._maxPtDist},Fi.prototype.interfaces_=function(){return[Ut]},Fi.prototype.getClass=function(){return Fi};var Gi=function(t,e,n){this._minValidDistance=null,this._maxValidDistance=null,this._minDistanceFound=null,this._maxDistanceFound=null,this._isValid=!0,this._errMsg=null,this._errorLocation=null,this._errorIndicator=null,this._input=t||null,this._bufDistance=e||null,this._result=n||null},qi={VERBOSE:{configurable:!0},MAX_DISTANCE_DIFF_FRAC:{configurable:!0}};Gi.prototype.checkMaximumDistance=function(t,e,n){var i=new Di(e,t);if(i.setDensifyFraction(.25),this._maxDistanceFound=i.orientedDistance(),this._maxDistanceFound>n){this._isValid=!1;var r=i.getCoordinates();this._errorLocation=r[1],this._errorIndicator=t.getFactory().createLineString(r),this._errMsg="Distance between buffer curve and input is too large ("+this._maxDistanceFound+" at "+Z.toLineString(r[0],r[1])+")"}},Gi.prototype.isValid=function(){var t=Math.abs(this._bufDistance),e=Gi.MAX_DISTANCE_DIFF_FRAC*t;return this._minValidDistance=t-e,this._maxValidDistance=t+e,!(!this._input.isEmpty()&&!this._result.isEmpty())||(this._bufDistance>0?this.checkPositiveValid():this.checkNegativeValid(),Gi.VERBOSE&&Y.out.println("Min Dist= "+this._minDistanceFound+"  err= "+(1-this._minDistanceFound/this._bufDistance)+"  Max Dist= "+this._maxDistanceFound+"  err= "+(this._maxDistanceFound/this._bufDistance-1)),this._isValid)},Gi.prototype.checkNegativeValid=function(){if(!(this._input instanceof $t||this._input instanceof ne||this._input instanceof zt))return null;var t=this.getPolygonLines(this._input);if(this.checkMinimumDistance(t,this._result,this._minValidDistance),!this._isValid)return null;this.checkMaximumDistance(t,this._result,this._maxValidDistance)},Gi.prototype.getErrorIndicator=function(){return this._errorIndicator},Gi.prototype.checkMinimumDistance=function(t,e,n){var i=new Ti(t,e,n);if(this._minDistanceFound=i.distance(),this._minDistanceFound<n){this._isValid=!1;var r=i.nearestPoints();this._errorLocation=i.nearestPoints()[1],this._errorIndicator=t.getFactory().createLineString(r),this._errMsg="Distance between buffer curve and input is too small ("+this._minDistanceFound+" at "+Z.toLineString(r[0],r[1])+" )"}},Gi.prototype.checkPositiveValid=function(){var t=this._result.getBoundary();if(this.checkMinimumDistance(this._input,t,this._minValidDistance),!this._isValid)return null;this.checkMaximumDistance(this._input,t,this._maxValidDistance)},Gi.prototype.getErrorLocation=function(){return this._errorLocation},Gi.prototype.getPolygonLines=function(t){for(var e=new Nt,n=new Ci(e),i=Ni.getPolygons(t).iterator();i.hasNext();){i.next().apply(n)}return t.getFactory().buildGeometry(e)},Gi.prototype.getErrorMessage=function(){return this._errMsg},Gi.prototype.interfaces_=function(){return[]},Gi.prototype.getClass=function(){return Gi},qi.VERBOSE.get=function(){return!1},qi.MAX_DISTANCE_DIFF_FRAC.get=function(){return.012},Object.defineProperties(Gi,qi);var Bi=function(t,e,n){this._isValid=!0,this._errorMsg=null,this._errorLocation=null,this._errorIndicator=null,this._input=t||null,this._distance=e||null,this._result=n||null},Vi={VERBOSE:{configurable:!0},MAX_ENV_DIFF_FRAC:{configurable:!0}};Bi.prototype.isValid=function(){return this.checkPolygonal(),this._isValid?(this.checkExpectedEmpty(),this._isValid?(this.checkEnvelope(),this._isValid?(this.checkArea(),this._isValid?(this.checkDistance(),this._isValid):this._isValid):this._isValid):this._isValid):this._isValid},Bi.prototype.checkEnvelope=function(){if(this._distance<0)return null;var t=this._distance*Bi.MAX_ENV_DIFF_FRAC;0===t&&(t=.001);var e=new j(this._input.getEnvelopeInternal());e.expandBy(this._distance);var n=new j(this._result.getEnvelopeInternal());n.expandBy(t),n.contains(e)||(this._isValid=!1,this._errorMsg="Buffer envelope is incorrect",this._errorIndicator=this._input.getFactory().toGeometry(n)),this.report("Envelope")},Bi.prototype.checkDistance=function(){var t=new Gi(this._input,this._distance,this._result);t.isValid()||(this._isValid=!1,this._errorMsg=t.getErrorMessage(),this._errorLocation=t.getErrorLocation(),this._errorIndicator=t.getErrorIndicator()),this.report("Distance")},Bi.prototype.checkArea=function(){var t=this._input.getArea(),e=this._result.getArea();this._distance>0&&t>e&&(this._isValid=!1,this._errorMsg="Area of positive buffer is smaller than input",this._errorIndicator=this._result),this._distance<0&&t<e&&(this._isValid=!1,this._errorMsg="Area of negative buffer is larger than input",this._errorIndicator=this._result),this.report("Area")},Bi.prototype.checkPolygonal=function(){this._result instanceof $t||this._result instanceof ne||(this._isValid=!1),this._errorMsg="Result is not polygonal",this._errorIndicator=this._result,this.report("Polygonal")},Bi.prototype.getErrorIndicator=function(){return this._errorIndicator},Bi.prototype.getErrorLocation=function(){return this._errorLocation},Bi.prototype.checkExpectedEmpty=function(){return this._input.getDimension()>=2?null:this._distance>0?null:(this._result.isEmpty()||(this._isValid=!1,this._errorMsg="Result is non-empty",this._errorIndicator=this._result),void this.report("ExpectedEmpty"))},Bi.prototype.report=function(t){if(!Bi.VERBOSE)return null;Y.out.println("Check "+t+": "+(this._isValid?"passed":"FAILED"))},Bi.prototype.getErrorMessage=function(){return this._errorMsg},Bi.prototype.interfaces_=function(){return[]},Bi.prototype.getClass=function(){return Bi},Bi.isValidMsg=function(t,e,n){var i=new Bi(t,e,n);return i.isValid()?null:i.getErrorMessage()},Bi.isValid=function(t,e,n){return!!new Bi(t,e,n).isValid()},Vi.VERBOSE.get=function(){return!1},Vi.MAX_ENV_DIFF_FRAC.get=function(){return.012},Object.defineProperties(Bi,Vi);var Ui=function(){this._pts=null,this._data=null;var t=arguments[0],e=arguments[1];this._pts=t,this._data=e};Ui.prototype.getCoordinates=function(){return this._pts},Ui.prototype.size=function(){return this._pts.length},Ui.prototype.getCoordinate=function(t){return this._pts[t]},Ui.prototype.isClosed=function(){return this._pts[0].equals(this._pts[this._pts.length-1])},Ui.prototype.getSegmentOctant=function(t){return t===this._pts.length-1?-1:pn.octant(this.getCoordinate(t),this.getCoordinate(t+1))},Ui.prototype.setData=function(t){this._data=t},Ui.prototype.getData=function(){return this._data},Ui.prototype.toString=function(){return Z.toLineString(new ue(this._pts))},Ui.prototype.interfaces_=function(){return[hn]},Ui.prototype.getClass=function(){return Ui};var zi=function(){this._findAllIntersections=!1,this._isCheckEndSegmentsOnly=!1,this._li=null,this._interiorIntersection=null,this._intSegments=null,this._intersections=new Nt,this._intersectionCount=0,this._keepIntersections=!0;var t=arguments[0];this._li=t,this._interiorIntersection=null};zi.prototype.getInteriorIntersection=function(){return this._interiorIntersection},zi.prototype.setCheckEndSegmentsOnly=function(t){this._isCheckEndSegmentsOnly=t},zi.prototype.getIntersectionSegments=function(){return this._intSegments},zi.prototype.count=function(){return this._intersectionCount},zi.prototype.getIntersections=function(){return this._intersections},zi.prototype.setFindAllIntersections=function(t){this._findAllIntersections=t},zi.prototype.setKeepIntersections=function(t){this._keepIntersections=t},zi.prototype.processIntersections=function(t,e,n,i){if(!this._findAllIntersections&&this.hasIntersection())return null;if(t===n&&e===i)return null;if(this._isCheckEndSegmentsOnly){if(!(this.isEndSegment(t,e)||this.isEndSegment(n,i)))return null}var r=t.getCoordinates()[e],o=t.getCoordinates()[e+1],s=n.getCoordinates()[i],a=n.getCoordinates()[i+1];this._li.computeIntersection(r,o,s,a),this._li.hasIntersection()&&this._li.isInteriorIntersection()&&(this._intSegments=new Array(4).fill(null),this._intSegments[0]=r,this._intSegments[1]=o,this._intSegments[2]=s,this._intSegments[3]=a,this._interiorIntersection=this._li.getIntersection(0),this._keepIntersections&&this._intersections.add(this._interiorIntersection),this._intersectionCount++)},zi.prototype.isEndSegment=function(t,e){return 0===e||e>=t.size()-2},zi.prototype.hasIntersection=function(){return null!==this._interiorIntersection},zi.prototype.isDone=function(){return!this._findAllIntersections&&null!==this._interiorIntersection},zi.prototype.interfaces_=function(){return[Wn]},zi.prototype.getClass=function(){return zi},zi.createAllIntersectionsFinder=function(t){var e=new zi(t);return e.setFindAllIntersections(!0),e},zi.createAnyIntersectionFinder=function(t){return new zi(t)},zi.createIntersectionCounter=function(t){var e=new zi(t);return e.setFindAllIntersections(!0),e.setKeepIntersections(!1),e};var Xi=function(){this._li=new rt,this._segStrings=null,this._findAllIntersections=!1,this._segInt=null,this._isValid=!0;var t=arguments[0];this._segStrings=t};Xi.prototype.execute=function(){if(null!==this._segInt)return null;this.checkInteriorIntersections()},Xi.prototype.getIntersections=function(){return this._segInt.getIntersections()},Xi.prototype.isValid=function(){return this.execute(),this._isValid},Xi.prototype.setFindAllIntersections=function(t){this._findAllIntersections=t},Xi.prototype.checkInteriorIntersections=function(){this._isValid=!0,this._segInt=new zi(this._li),this._segInt.setFindAllIntersections(this._findAllIntersections);var t=new xn;if(t.setSegmentIntersector(this._segInt),t.computeNodes(this._segStrings),this._segInt.hasIntersection())return this._isValid=!1,null},Xi.prototype.checkValid=function(){if(this.execute(),!this._isValid)throw new we(this.getErrorMessage(),this._segInt.getInteriorIntersection())},Xi.prototype.getErrorMessage=function(){if(this._isValid)return"no intersections found";var t=this._segInt.getIntersectionSegments();return"found non-noded intersection between "+Z.toLineString(t[0],t[1])+" and "+Z.toLineString(t[2],t[3])},Xi.prototype.interfaces_=function(){return[]},Xi.prototype.getClass=function(){return Xi},Xi.computeIntersections=function(t){var e=new Xi(t);return e.setFindAllIntersections(!0),e.isValid(),e.getIntersections()};var Yi=function t(){this._nv=null;var e=arguments[0];this._nv=new Xi(t.toSegmentStrings(e))};Yi.prototype.checkValid=function(){this._nv.checkValid()},Yi.prototype.interfaces_=function(){return[]},Yi.prototype.getClass=function(){return Yi},Yi.toSegmentStrings=function(t){for(var e=new Nt,n=t.iterator();n.hasNext();){var i=n.next();e.add(new Ui(i.getCoordinates(),i))}return e},Yi.checkValid=function(t){new Yi(t).checkValid()};var ki=function(t){this._mapOp=t};ki.prototype.map=function(t){for(var e=new Nt,n=0;n<t.getNumGeometries();n++){var i=this._mapOp.map(t.getGeometryN(n));i.isEmpty()||e.add(i)}return t.getFactory().createGeometryCollection(_e.toGeometryArray(e))},ki.prototype.interfaces_=function(){return[]},ki.prototype.getClass=function(){return ki},ki.map=function(t,e){return new ki(e).map(t)};var ji=function(){this._op=null,this._geometryFactory=null,this._ptLocator=null,this._lineEdgesList=new Nt,this._resultLineList=new Nt;var t=arguments[0],e=arguments[1],n=arguments[2];this._op=t,this._geometryFactory=e,this._ptLocator=n};ji.prototype.collectLines=function(t){for(var e=this._op.getGraph().getEdgeEnds().iterator();e.hasNext();){var n=e.next();this.collectLineEdge(n,t,this._lineEdgesList),this.collectBoundaryTouchEdge(n,t,this._lineEdgesList)}},ji.prototype.labelIsolatedLine=function(t,e){var n=this._ptLocator.locate(t.getCoordinate(),this._op.getArgGeometry(e));t.getLabel().setLocation(e,n)},ji.prototype.build=function(t){return this.findCoveredLineEdges(),this.collectLines(t),this.buildLines(t),this._resultLineList},ji.prototype.collectLineEdge=function(t,e,n){var i=t.getLabel(),r=t.getEdge();t.isLineEdge()&&(t.isVisited()||!Lr.isResultOfOp(i,e)||r.isCovered()||(n.add(r),t.setVisitedEdge(!0)))},ji.prototype.findCoveredLineEdges=function(){for(var t=this._op.getGraph().getNodes().iterator();t.hasNext();){t.next().getEdges().findCoveredLineEdges()}for(var e=this._op.getGraph().getEdgeEnds().iterator();e.hasNext();){var n=e.next(),i=n.getEdge();if(n.isLineEdge()&&!i.isCoveredSet()){var r=this._op.isCoveredByA(n.getCoordinate());i.setCovered(r)}}},ji.prototype.labelIsolatedLines=function(t){for(var e=t.iterator();e.hasNext();){var n=e.next(),i=n.getLabel();n.isIsolated()&&(i.isNull(0)?this.labelIsolatedLine(n,0):this.labelIsolatedLine(n,1))}},ji.prototype.buildLines=function(t){for(var e=this._lineEdgesList.iterator();e.hasNext();){var n=e.next(),i=this._geometryFactory.createLineString(n.getCoordinates());this._resultLineList.add(i),n.setInResult(!0)}},ji.prototype.collectBoundaryTouchEdge=function(t,e,n){var i=t.getLabel();return t.isLineEdge()?null:t.isVisited()?null:t.isInteriorAreaEdge()?null:t.getEdge().isInResult()?null:(et.isTrue(!(t.isInResult()||t.getSym().isInResult())||!t.getEdge().isInResult()),void(Lr.isResultOfOp(i,e)&&e===Lr.INTERSECTION&&(n.add(t.getEdge()),t.setVisitedEdge(!0))))},ji.prototype.interfaces_=function(){return[]},ji.prototype.getClass=function(){return ji};var Hi=function(){this._op=null,this._geometryFactory=null,this._resultPointList=new Nt;var t=arguments[0],e=arguments[1];this._op=t,this._geometryFactory=e};Hi.prototype.filterCoveredNodeToPoint=function(t){var e=t.getCoordinate();if(!this._op.isCoveredByLA(e)){var n=this._geometryFactory.createPoint(e);this._resultPointList.add(n)}},Hi.prototype.extractNonCoveredResultNodes=function(t){for(var e=this._op.getGraph().getNodes().iterator();e.hasNext();){var n=e.next();if(!n.isInResult()&&(!n.isIncidentEdgeInResult()&&(0===n.getEdges().getDegree()||t===Lr.INTERSECTION))){var i=n.getLabel();Lr.isResultOfOp(i,t)&&this.filterCoveredNodeToPoint(n)}}},Hi.prototype.build=function(t){return this.extractNonCoveredResultNodes(t),this._resultPointList},Hi.prototype.interfaces_=function(){return[]},Hi.prototype.getClass=function(){return Hi};var Wi=function(){this._inputGeom=null,this._factory=null,this._pruneEmptyGeometry=!0,this._preserveGeometryCollectionType=!0,this._preserveCollections=!1,this._preserveType=!1};Wi.prototype.transformPoint=function(t,e){return this._factory.createPoint(this.transformCoordinates(t.getCoordinateSequence(),t))},Wi.prototype.transformPolygon=function(t,e){var n=!0,i=this.transformLinearRing(t.getExteriorRing(),t);null!==i&&i instanceof ee&&!i.isEmpty()||(n=!1);for(var r=new Nt,o=0;o<t.getNumInteriorRing();o++){var s=this.transformLinearRing(t.getInteriorRingN(o),t);null===s||s.isEmpty()||(s instanceof ee||(n=!1),r.add(s))}if(n)return this._factory.createPolygon(i,r.toArray([]));var a=new Nt;return null!==i&&a.add(i),a.addAll(r),this._factory.buildGeometry(a)},Wi.prototype.createCoordinateSequence=function(t){return this._factory.getCoordinateSequenceFactory().create(t)},Wi.prototype.getInputGeometry=function(){return this._inputGeom},Wi.prototype.transformMultiLineString=function(t,e){for(var n=new Nt,i=0;i<t.getNumGeometries();i++){var r=this.transformLineString(t.getGeometryN(i),t);null!==r&&(r.isEmpty()||n.add(r))}return this._factory.buildGeometry(n)},Wi.prototype.transformCoordinates=function(t,e){return this.copy(t)},Wi.prototype.transformLineString=function(t,e){return this._factory.createLineString(this.transformCoordinates(t.getCoordinateSequence(),t))},Wi.prototype.transformMultiPoint=function(t,e){for(var n=new Nt,i=0;i<t.getNumGeometries();i++){var r=this.transformPoint(t.getGeometryN(i),t);null!==r&&(r.isEmpty()||n.add(r))}return this._factory.buildGeometry(n)},Wi.prototype.transformMultiPolygon=function(t,e){for(var n=new Nt,i=0;i<t.getNumGeometries();i++){var r=this.transformPolygon(t.getGeometryN(i),t);null!==r&&(r.isEmpty()||n.add(r))}return this._factory.buildGeometry(n)},Wi.prototype.copy=function(t){return t.copy()},Wi.prototype.transformGeometryCollection=function(t,e){for(var n=new Nt,i=0;i<t.getNumGeometries();i++){var r=this.transform(t.getGeometryN(i));null!==r&&(this._pruneEmptyGeometry&&r.isEmpty()||n.add(r))}return this._preserveGeometryCollectionType?this._factory.createGeometryCollection(_e.toGeometryArray(n)):this._factory.buildGeometry(n)},Wi.prototype.transform=function(t){if(this._inputGeom=t,this._factory=t.getFactory(),t instanceof Qt)return this.transformPoint(t,null);if(t instanceof te)return this.transformMultiPoint(t,null);if(t instanceof ee)return this.transformLinearRing(t,null);if(t instanceof Kt)return this.transformLineString(t,null);if(t instanceof Xt)return this.transformMultiLineString(t,null);if(t instanceof $t)return this.transformPolygon(t,null);if(t instanceof ne)return this.transformMultiPolygon(t,null);if(t instanceof zt)return this.transformGeometryCollection(t,null);throw new m("Unknown Geometry subtype: "+t.getClass().getName())},Wi.prototype.transformLinearRing=function(t,e){var n=this.transformCoordinates(t.getCoordinateSequence(),t);if(null===n)return this._factory.createLinearRing(null);var i=n.size();return i>0&&i<4&&!this._preserveType?this._factory.createLineString(n):this._factory.createLinearRing(n)},Wi.prototype.interfaces_=function(){return[]},Wi.prototype.getClass=function(){return Wi};var Ki=function t(){if(this._snapTolerance=0,this._srcPts=null,this._seg=new dn,this._allowSnappingToSourceVertices=!1,this._isClosed=!1,arguments[0]instanceof Kt&&"number"==typeof arguments[1]){var e=arguments[0],n=arguments[1];t.call(this,e.getCoordinates(),n)}else if(arguments[0]instanceof Array&&"number"==typeof arguments[1]){var i=arguments[0],r=arguments[1];this._srcPts=i,this._isClosed=t.isClosed(i),this._snapTolerance=r}};Ki.prototype.snapVertices=function(t,e){for(var n=this._isClosed?t.size()-1:t.size(),i=0;i<n;i++){var r=t.get(i),o=this.findSnapForVertex(r,e);null!==o&&(t.set(i,new C(o)),0===i&&this._isClosed&&t.set(t.size()-1,new C(o)))}},Ki.prototype.findSnapForVertex=function(t,e){for(var n=0;n<e.length;n++){if(t.equals2D(e[n]))return null;if(t.distance(e[n])<this._snapTolerance)return e[n]}return null},Ki.prototype.snapTo=function(t){var e=new St(this._srcPts);this.snapVertices(e,t),this.snapSegments(e,t);return e.toCoordinateArray()},Ki.prototype.snapSegments=function(t,e){if(0===e.length)return null;var n=e.length;e[0].equals2D(e[e.length-1])&&(n=e.length-1);for(var i=0;i<n;i++){var r=e[i],o=this.findSegmentIndexToSnap(r,t);o>=0&&t.add(o+1,new C(r),!1)}},Ki.prototype.findSegmentIndexToSnap=function(t,e){for(var n=v.MAX_VALUE,i=-1,r=0;r<e.size()-1;r++){if(this._seg.p0=e.get(r),this._seg.p1=e.get(r+1),this._seg.p0.equals2D(t)||this._seg.p1.equals2D(t)){if(this._allowSnappingToSourceVertices)continue;return-1}var o=this._seg.distance(t);o<this._snapTolerance&&o<n&&(n=o,i=r)}return i},Ki.prototype.setAllowSnappingToSourceVertices=function(t){this._allowSnappingToSourceVertices=t},Ki.prototype.interfaces_=function(){return[]},Ki.prototype.getClass=function(){return Ki},Ki.isClosed=function(t){return!(t.length<=1)&&t[0].equals2D(t[t.length-1])};var Ji=function(t){this._srcGeom=t||null},Qi={SNAP_PRECISION_FACTOR:{configurable:!0}};Ji.prototype.snapTo=function(t,e){var n=this.extractTargetCoordinates(t);return new Zi(e,n).transform(this._srcGeom)},Ji.prototype.snapToSelf=function(t,e){var n=this.extractTargetCoordinates(this._srcGeom),i=new Zi(t,n,!0).transform(this._srcGeom),r=i;return e&&T(r,Zt)&&(r=i.buffer(0)),r},Ji.prototype.computeSnapTolerance=function(t){return this.computeMinimumSegmentLength(t)/10},Ji.prototype.extractTargetCoordinates=function(t){for(var e=new f,n=t.getCoordinates(),i=0;i<n.length;i++)e.add(n[i]);return e.toArray(new Array(0).fill(null))},Ji.prototype.computeMinimumSegmentLength=function(t){for(var e=v.MAX_VALUE,n=0;n<t.length-1;n++){var i=t[n].distance(t[n+1]);i<e&&(e=i)}return e},Ji.prototype.interfaces_=function(){return[]},Ji.prototype.getClass=function(){return Ji},Ji.snap=function(t,e,n){var i=new Array(2).fill(null),r=new Ji(t);i[0]=r.snapTo(e,n);var o=new Ji(e);return i[1]=o.snapTo(i[0],n),i},Ji.computeOverlaySnapTolerance=function(){if(1===arguments.length){var t=arguments[0],e=Ji.computeSizeBasedSnapTolerance(t),n=t.getPrecisionModel();if(n.getType()===fe.FIXED){var i=1/n.getScale()*2/1.415;i>e&&(e=i)}return e}if(2===arguments.length){var r=arguments[0],o=arguments[1];return Math.min(Ji.computeOverlaySnapTolerance(r),Ji.computeOverlaySnapTolerance(o))}},Ji.computeSizeBasedSnapTolerance=function(t){var e=t.getEnvelopeInternal();return Math.min(e.getHeight(),e.getWidth())*Ji.SNAP_PRECISION_FACTOR},Ji.snapToSelf=function(t,e,n){return new Ji(t).snapToSelf(e,n)},Qi.SNAP_PRECISION_FACTOR.get=function(){return 1e-9},Object.defineProperties(Ji,Qi);var Zi=function(t){function e(e,n,i){t.call(this),this._snapTolerance=e||null,this._snapPts=n||null,this._isSelfSnap=void 0!==i&&i}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.snapLine=function(t,e){var n=new Ki(t,this._snapTolerance);return n.setAllowSnappingToSourceVertices(this._isSelfSnap),n.snapTo(e)},e.prototype.transformCoordinates=function(t,e){var n=t.toCoordinateArray(),i=this.snapLine(n,this._snapPts);return this._factory.getCoordinateSequenceFactory().create(i)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Wi),$i=function(){this._isFirst=!0,this._commonMantissaBitsCount=53,this._commonBits=0,this._commonSignExp=null};$i.prototype.getCommon=function(){return v.longBitsToDouble(this._commonBits)},$i.prototype.add=function(t){var e=v.doubleToLongBits(t);if(this._isFirst)return this._commonBits=e,this._commonSignExp=$i.signExpBits(this._commonBits),this._isFirst=!1,null;if($i.signExpBits(e)!==this._commonSignExp)return this._commonBits=0,null;this._commonMantissaBitsCount=$i.numCommonMostSigMantissaBits(this._commonBits,e),this._commonBits=$i.zeroLowerBits(this._commonBits,64-(12+this._commonMantissaBitsCount))},$i.prototype.toString=function(){if(1===arguments.length){var t=arguments[0],e=v.longBitsToDouble(t),n="0000000000000000000000000000000000000000000000000000000000000000"+v.toBinaryString(t),i=n.substring(n.length-64);return i.substring(0,1)+"  "+i.substring(1,12)+"(exp) "+i.substring(12)+" [ "+e+" ]"}},$i.prototype.interfaces_=function(){return[]},$i.prototype.getClass=function(){return $i},$i.getBit=function(t,e){return 0!=(t&1<<e)?1:0},$i.signExpBits=function(t){return t>>52},$i.zeroLowerBits=function(t,e){return t&~((1<<e)-1)},$i.numCommonMostSigMantissaBits=function(t,e){for(var n=0,i=52;i>=0;i--){if($i.getBit(t,i)!==$i.getBit(e,i))return n;n++}return 52};var tr=function(){this._commonCoord=null,this._ccFilter=new nr},er={CommonCoordinateFilter:{configurable:!0},Translater:{configurable:!0}};tr.prototype.addCommonBits=function(t){var e=new ir(this._commonCoord);t.apply(e),t.geometryChanged()},tr.prototype.removeCommonBits=function(t){if(0===this._commonCoord.x&&0===this._commonCoord.y)return t;var e=new C(this._commonCoord);e.x=-e.x,e.y=-e.y;var n=new ir(e);return t.apply(n),t.geometryChanged(),t},tr.prototype.getCommonCoordinate=function(){return this._commonCoord},tr.prototype.add=function(t){t.apply(this._ccFilter),this._commonCoord=this._ccFilter.getCommonCoordinate()},tr.prototype.interfaces_=function(){return[]},tr.prototype.getClass=function(){return tr},er.CommonCoordinateFilter.get=function(){return nr},er.Translater.get=function(){return ir},Object.defineProperties(tr,er);var nr=function(){this._commonBitsX=new $i,this._commonBitsY=new $i};nr.prototype.filter=function(t){this._commonBitsX.add(t.x),this._commonBitsY.add(t.y)},nr.prototype.getCommonCoordinate=function(){return new C(this._commonBitsX.getCommon(),this._commonBitsY.getCommon())},nr.prototype.interfaces_=function(){return[ft]},nr.prototype.getClass=function(){return nr};var ir=function(){this.trans=null;var t=arguments[0];this.trans=t};ir.prototype.filter=function(t,e){var n=t.getOrdinate(e,0)+this.trans.x,i=t.getOrdinate(e,1)+this.trans.y;t.setOrdinate(e,0,n),t.setOrdinate(e,1,i)},ir.prototype.isDone=function(){return!1},ir.prototype.isGeometryChanged=function(){return!0},ir.prototype.interfaces_=function(){return[Ut]},ir.prototype.getClass=function(){return ir};var rr=function(t,e){this._geom=new Array(2).fill(null),this._snapTolerance=null,this._cbr=null,this._geom[0]=t,this._geom[1]=e,this.computeSnapTolerance()};rr.prototype.selfSnap=function(t){return new Ji(t).snapTo(t,this._snapTolerance)},rr.prototype.removeCommonBits=function(t){this._cbr=new tr,this._cbr.add(t[0]),this._cbr.add(t[1]);var e=new Array(2).fill(null);return e[0]=this._cbr.removeCommonBits(t[0].copy()),e[1]=this._cbr.removeCommonBits(t[1].copy()),e},rr.prototype.prepareResult=function(t){return this._cbr.addCommonBits(t),t},rr.prototype.getResultGeometry=function(t){var e=this.snap(this._geom),n=Lr.overlayOp(e[0],e[1],t);return this.prepareResult(n)},rr.prototype.checkValid=function(t){t.isValid()||Y.out.println("Snapped geometry is invalid")},rr.prototype.computeSnapTolerance=function(){this._snapTolerance=Ji.computeOverlaySnapTolerance(this._geom[0],this._geom[1])},rr.prototype.snap=function(t){var e=this.removeCommonBits(t);return Ji.snap(e[0],e[1],this._snapTolerance)},rr.prototype.interfaces_=function(){return[]},rr.prototype.getClass=function(){return rr},rr.overlayOp=function(t,e,n){return new rr(t,e).getResultGeometry(n)},rr.union=function(t,e){return rr.overlayOp(t,e,Lr.UNION)},rr.intersection=function(t,e){return rr.overlayOp(t,e,Lr.INTERSECTION)},rr.symDifference=function(t,e){return rr.overlayOp(t,e,Lr.SYMDIFFERENCE)},rr.difference=function(t,e){return rr.overlayOp(t,e,Lr.DIFFERENCE)};var or=function(t,e){this._geom=new Array(2).fill(null),this._geom[0]=t,this._geom[1]=e};or.prototype.getResultGeometry=function(t){var e=null,n=!1,i=null;try{e=Lr.overlayOp(this._geom[0],this._geom[1],t);n=!0}catch(t){if(!(t instanceof $))throw t;i=t}if(!n)try{e=rr.overlayOp(this._geom[0],this._geom[1],t)}catch(t){throw t instanceof $?i:t}return e},or.prototype.interfaces_=function(){return[]},or.prototype.getClass=function(){return or},or.overlayOp=function(t,e,n){return new or(t,e).getResultGeometry(n)},or.union=function(t,e){return or.overlayOp(t,e,Lr.UNION)},or.intersection=function(t,e){return or.overlayOp(t,e,Lr.INTERSECTION)},or.symDifference=function(t,e){return or.overlayOp(t,e,Lr.SYMDIFFERENCE)},or.difference=function(t,e){return or.overlayOp(t,e,Lr.DIFFERENCE)};var sr=function(){this.mce=null,this.chainIndex=null;var t=arguments[0],e=arguments[1];this.mce=t,this.chainIndex=e};sr.prototype.computeIntersections=function(t,e){this.mce.computeIntersectsForChain(this.chainIndex,t.mce,t.chainIndex,e)},sr.prototype.interfaces_=function(){return[]},sr.prototype.getClass=function(){return sr};var ar=function t(){if(this._label=null,this._xValue=null,this._eventType=null,this._insertEvent=null,this._deleteEventIndex=null,this._obj=null,2===arguments.length){var e=arguments[0],n=arguments[1];this._eventType=t.DELETE,this._xValue=e,this._insertEvent=n}else if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2];this._eventType=t.INSERT,this._label=i,this._xValue=r,this._obj=o}},ur={INSERT:{configurable:!0},DELETE:{configurable:!0}};ar.prototype.isDelete=function(){return this._eventType===ar.DELETE},ar.prototype.setDeleteEventIndex=function(t){this._deleteEventIndex=t},ar.prototype.getObject=function(){return this._obj},ar.prototype.compareTo=function(t){var e=t;return this._xValue<e._xValue?-1:this._xValue>e._xValue?1:this._eventType<e._eventType?-1:this._eventType>e._eventType?1:0},ar.prototype.getInsertEvent=function(){return this._insertEvent},ar.prototype.isInsert=function(){return this._eventType===ar.INSERT},ar.prototype.isSameLabel=function(t){return null!==this._label&&this._label===t._label},ar.prototype.getDeleteEventIndex=function(){return this._deleteEventIndex},ar.prototype.interfaces_=function(){return[E]},ar.prototype.getClass=function(){return ar},ur.INSERT.get=function(){return 1},ur.DELETE.get=function(){return 2},Object.defineProperties(ar,ur);var lr=function(){};lr.prototype.interfaces_=function(){return[]},lr.prototype.getClass=function(){return lr};var cr=function(){this._hasIntersection=!1,this._hasProper=!1,this._hasProperInterior=!1,this._properIntersectionPoint=null,this._li=null,this._includeProper=null,this._recordIsolated=null,this._isSelfIntersection=null,this._numIntersections=0,this.numTests=0,this._bdyNodes=null,this._isDone=!1,this._isDoneWhenProperInt=!1;var t=arguments[0],e=arguments[1],n=arguments[2];this._li=t,this._includeProper=e,this._recordIsolated=n};cr.prototype.isTrivialIntersection=function(t,e,n,i){if(t===n&&1===this._li.getIntersectionNum()){if(cr.isAdjacentSegments(e,i))return!0;if(t.isClosed()){var r=t.getNumPoints()-1;if(0===e&&i===r||0===i&&e===r)return!0}}return!1},cr.prototype.getProperIntersectionPoint=function(){return this._properIntersectionPoint},cr.prototype.setIsDoneIfProperInt=function(t){this._isDoneWhenProperInt=t},cr.prototype.hasProperInteriorIntersection=function(){return this._hasProperInterior},cr.prototype.isBoundaryPointInternal=function(t,e){for(var n=e.iterator();n.hasNext();){var i=n.next().getCoordinate();if(t.isIntersection(i))return!0}return!1},cr.prototype.hasProperIntersection=function(){return this._hasProper},cr.prototype.hasIntersection=function(){return this._hasIntersection},cr.prototype.isDone=function(){return this._isDone},cr.prototype.isBoundaryPoint=function(t,e){return null!==e&&(!!this.isBoundaryPointInternal(t,e[0])||!!this.isBoundaryPointInternal(t,e[1]))},cr.prototype.setBoundaryNodes=function(t,e){this._bdyNodes=new Array(2).fill(null),this._bdyNodes[0]=t,this._bdyNodes[1]=e},cr.prototype.addIntersections=function(t,e,n,i){if(t===n&&e===i)return null;this.numTests++;var r=t.getCoordinates()[e],o=t.getCoordinates()[e+1],s=n.getCoordinates()[i],a=n.getCoordinates()[i+1];this._li.computeIntersection(r,o,s,a),this._li.hasIntersection()&&(this._recordIsolated&&(t.setIsolated(!1),n.setIsolated(!1)),this._numIntersections++,this.isTrivialIntersection(t,e,n,i)||(this._hasIntersection=!0,!this._includeProper&&this._li.isProper()||(t.addIntersections(this._li,e,0),n.addIntersections(this._li,i,1)),this._li.isProper()&&(this._properIntersectionPoint=this._li.getIntersection(0).copy(),this._hasProper=!0,this._isDoneWhenProperInt&&(this._isDone=!0),this.isBoundaryPoint(this._li,this._bdyNodes)||(this._hasProperInterior=!0))))},cr.prototype.interfaces_=function(){return[]},cr.prototype.getClass=function(){return cr},cr.isAdjacentSegments=function(t,e){return 1===Math.abs(t-e)};var pr=function(t){function e(){t.call(this),this.events=new Nt,this.nOverlaps=null}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.prepareEvents=function(){$e.sort(this.events);for(var t=0;t<this.events.size();t++){var e=this.events.get(t);e.isDelete()&&e.getInsertEvent().setDeleteEventIndex(t)}},e.prototype.computeIntersections=function(){if(1===arguments.length){var t=arguments[0];this.nOverlaps=0,this.prepareEvents();for(var e=0;e<this.events.size();e++){var n=this.events.get(e);if(n.isInsert()&&this.processOverlaps(e,n.getDeleteEventIndex(),n,t),t.isDone())break}}else if(3===arguments.length)if(arguments[2]instanceof cr&&T(arguments[0],xt)&&T(arguments[1],xt)){var i=arguments[0],r=arguments[1],o=arguments[2];this.addEdges(i,i),this.addEdges(r,r),this.computeIntersections(o)}else if("boolean"==typeof arguments[2]&&T(arguments[0],xt)&&arguments[1]instanceof cr){var s=arguments[0],a=arguments[1];arguments[2]?this.addEdges(s,null):this.addEdges(s),this.computeIntersections(a)}},e.prototype.addEdge=function(t,e){for(var n=t.getMonotoneChainEdge(),i=n.getStartIndexes(),r=0;r<i.length-1;r++){var o=new sr(n,r),s=new ar(e,n.getMinX(r),o);this.events.add(s),this.events.add(new ar(n.getMaxX(r),s))}},e.prototype.processOverlaps=function(t,e,n,i){for(var r=n.getObject(),o=t;o<e;o++){var s=this.events.get(o);if(s.isInsert()){var a=s.getObject();n.isSameLabel(s)||(r.computeIntersections(a,i),this.nOverlaps++)}}},e.prototype.addEdges=function(){if(1===arguments.length)for(var t=arguments[0].iterator();t.hasNext();){var e=t.next();this.addEdge(e,e)}else if(2===arguments.length)for(var n=arguments[0],i=arguments[1],r=n.iterator();r.hasNext();){var o=r.next();this.addEdge(o,i)}},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(lr),hr=function(){this._min=v.POSITIVE_INFINITY,this._max=v.NEGATIVE_INFINITY},fr={NodeComparator:{configurable:!0}};hr.prototype.getMin=function(){return this._min},hr.prototype.intersects=function(t,e){return!(this._min>e||this._max<t)},hr.prototype.getMax=function(){return this._max},hr.prototype.toString=function(){return Z.toLineString(new C(this._min,0),new C(this._max,0))},hr.prototype.interfaces_=function(){return[]},hr.prototype.getClass=function(){return hr},fr.NodeComparator.get=function(){return gr},Object.defineProperties(hr,fr);var gr=function(){};gr.prototype.compare=function(t,e){var n=t,i=e,r=(n._min+n._max)/2,o=(i._min+i._max)/2;return r<o?-1:r>o?1:0},gr.prototype.interfaces_=function(){return[N]},gr.prototype.getClass=function(){return gr};var dr=function(t){function e(){t.call(this),this._item=null;var e=arguments[0],n=arguments[1],i=arguments[2];this._min=e,this._max=n,this._item=i}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.query=function(t,e,n){if(!this.intersects(t,e))return null;n.visitItem(this._item)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(hr),yr=function(t){function e(){t.call(this),this._node1=null,this._node2=null;var e=arguments[0],n=arguments[1];this._node1=e,this._node2=n,this.buildExtent(this._node1,this._node2)}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.buildExtent=function(t,e){this._min=Math.min(t._min,e._min),this._max=Math.max(t._max,e._max)},e.prototype.query=function(t,e,n){if(!this.intersects(t,e))return null;null!==this._node1&&this._node1.query(t,e,n),null!==this._node2&&this._node2.query(t,e,n)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(hr),_r=function(){this._leaves=new Nt,this._root=null,this._level=0};_r.prototype.buildTree=function(){$e.sort(this._leaves,new hr.NodeComparator);for(var t=this._leaves,e=null,n=new Nt;;){if(this.buildLevel(t,n),1===n.size())return n.get(0);e=t,t=n,n=e}},_r.prototype.insert=function(t,e,n){if(null!==this._root)throw new Error("Index cannot be added to once it has been queried");this._leaves.add(new dr(t,e,n))},_r.prototype.query=function(t,e,n){this.init(),this._root.query(t,e,n)},_r.prototype.buildRoot=function(){if(null!==this._root)return null;this._root=this.buildTree()},_r.prototype.printNode=function(t){Y.out.println(Z.toLineString(new C(t._min,this._level),new C(t._max,this._level)))},_r.prototype.init=function(){if(null!==this._root)return null;this.buildRoot()},_r.prototype.buildLevel=function(t,e){this._level++,e.clear();for(var n=0;n<t.size();n+=2){var i=t.get(n);if(null===(n+1<t.size()?t.get(n):null))e.add(i);else{var r=new yr(t.get(n),t.get(n+1));e.add(r)}}},_r.prototype.interfaces_=function(){return[]},_r.prototype.getClass=function(){return _r};var mr=function(){this._items=new Nt};mr.prototype.visitItem=function(t){this._items.add(t)},mr.prototype.getItems=function(){return this._items},mr.prototype.interfaces_=function(){return[Ke]},mr.prototype.getClass=function(){return mr};var vr=function(){this._index=null;var t=arguments[0];if(!T(t,Zt))throw new m("Argument must be Polygonal");this._index=new xr(t)},Ir={SegmentVisitor:{configurable:!0},IntervalIndexedGeometry:{configurable:!0}};vr.prototype.locate=function(t){var e=new st(t),n=new Er(e);return this._index.query(t.y,t.y,n),e.getLocation()},vr.prototype.interfaces_=function(){return[Vn]},vr.prototype.getClass=function(){return vr},Ir.SegmentVisitor.get=function(){return Er},Ir.IntervalIndexedGeometry.get=function(){return xr},Object.defineProperties(vr,Ir);var Er=function(){this._counter=null;var t=arguments[0];this._counter=t};Er.prototype.visitItem=function(t){var e=t;this._counter.countSegment(e.getCoordinate(0),e.getCoordinate(1))},Er.prototype.interfaces_=function(){return[Ke]},Er.prototype.getClass=function(){return Er};var xr=function(){this._index=new _r;var t=arguments[0];this.init(t)};xr.prototype.init=function(t){for(var e=Ci.getLines(t).iterator();e.hasNext();){var n=e.next().getCoordinates();this.addLine(n)}},xr.prototype.addLine=function(t){for(var e=1;e<t.length;e++){var n=new dn(t[e-1],t[e]),i=Math.min(n.p0.y,n.p1.y),r=Math.max(n.p0.y,n.p1.y);this._index.insert(i,r,n)}},xr.prototype.query=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1],n=new mr;return this._index.query(t,e,n),n.getItems()}if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2];this._index.query(i,r,o)}},xr.prototype.interfaces_=function(){return[]},xr.prototype.getClass=function(){return xr};var Nr=function(t){function e(){if(t.call(this),this._parentGeom=null,this._lineEdgeMap=new he,this._boundaryNodeRule=null,this._useBoundaryDeterminationRule=!0,this._argIndex=null,this._boundaryNodes=null,this._hasTooFewPoints=!1,this._invalidPoint=null,this._areaPtLocator=null,this._ptLocator=new Si,2===arguments.length){var e=arguments[0],n=arguments[1],i=gt.OGC_SFS_BOUNDARY_RULE;this._argIndex=e,this._parentGeom=n,this._boundaryNodeRule=i,null!==n&&this.add(n)}else if(3===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2];this._argIndex=r,this._parentGeom=o,this._boundaryNodeRule=s,null!==o&&this.add(o)}}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.insertBoundaryPoint=function(t,n){var i=this._nodes.addNode(n).getLabel(),r=1;w.NONE;i.getLocation(t,Se.ON)===w.BOUNDARY&&r++;var o=e.determineBoundary(this._boundaryNodeRule,r);i.setLocation(t,o)},e.prototype.computeSelfNodes=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1];return this.computeSelfNodes(t,e,!1)}if(3===arguments.length){var n=arguments[0],i=arguments[1],r=arguments[2],o=new cr(n,!0,!1);o.setIsDoneIfProperInt(r);var s=this.createEdgeSetIntersector(),a=this._parentGeom instanceof ee||this._parentGeom instanceof $t||this._parentGeom instanceof ne,u=i||!a;return s.computeIntersections(this._edges,o,u),this.addSelfIntersectionNodes(this._argIndex),o}},e.prototype.computeSplitEdges=function(t){for(var e=this._edges.iterator();e.hasNext();){e.next().eiList.addSplitEdges(t)}},e.prototype.computeEdgeIntersections=function(t,e,n){var i=new cr(e,n,!0);i.setBoundaryNodes(this.getBoundaryNodes(),t.getBoundaryNodes());return this.createEdgeSetIntersector().computeIntersections(this._edges,t._edges,i),i},e.prototype.getGeometry=function(){return this._parentGeom},e.prototype.getBoundaryNodeRule=function(){return this._boundaryNodeRule},e.prototype.hasTooFewPoints=function(){return this._hasTooFewPoints},e.prototype.addPoint=function(){if(arguments[0]instanceof Qt){var t=arguments[0].getCoordinate();this.insertPoint(this._argIndex,t,w.INTERIOR)}else if(arguments[0]instanceof C){var e=arguments[0];this.insertPoint(this._argIndex,e,w.INTERIOR)}},e.prototype.addPolygon=function(t){this.addPolygonRing(t.getExteriorRing(),w.EXTERIOR,w.INTERIOR);for(var e=0;e<t.getNumInteriorRing();e++){var n=t.getInteriorRingN(e);this.addPolygonRing(n,w.INTERIOR,w.EXTERIOR)}},e.prototype.addEdge=function(t){this.insertEdge(t);var e=t.getCoordinates();this.insertPoint(this._argIndex,e[0],w.BOUNDARY),this.insertPoint(this._argIndex,e[e.length-1],w.BOUNDARY)},e.prototype.addLineString=function(t){var e=Lt.removeRepeatedPoints(t.getCoordinates());if(e.length<2)return this._hasTooFewPoints=!0,this._invalidPoint=e[0],null;var n=new ni(e,new Pe(this._argIndex,w.INTERIOR));this._lineEdgeMap.put(t,n),this.insertEdge(n),et.isTrue(e.length>=2,"found LineString with single point"),this.insertBoundaryPoint(this._argIndex,e[0]),this.insertBoundaryPoint(this._argIndex,e[e.length-1])},e.prototype.getInvalidPoint=function(){return this._invalidPoint},e.prototype.getBoundaryPoints=function(){for(var t=this.getBoundaryNodes(),e=new Array(t.size()).fill(null),n=0,i=t.iterator();i.hasNext();){var r=i.next();e[n++]=r.getCoordinate().copy()}return e},e.prototype.getBoundaryNodes=function(){return null===this._boundaryNodes&&(this._boundaryNodes=this._nodes.getBoundaryNodes(this._argIndex)),this._boundaryNodes},e.prototype.addSelfIntersectionNode=function(t,e,n){if(this.isBoundaryNode(t,e))return null;n===w.BOUNDARY&&this._useBoundaryDeterminationRule?this.insertBoundaryPoint(t,e):this.insertPoint(t,e,n)},e.prototype.addPolygonRing=function(t,e,n){if(t.isEmpty())return null;var i=Lt.removeRepeatedPoints(t.getCoordinates());if(i.length<4)return this._hasTooFewPoints=!0,this._invalidPoint=i[0],null;var r=e,o=n;at.isCCW(i)&&(r=n,o=e);var s=new ni(i,new Pe(this._argIndex,w.BOUNDARY,r,o));this._lineEdgeMap.put(t,s),this.insertEdge(s),this.insertPoint(this._argIndex,i[0],w.BOUNDARY)},e.prototype.insertPoint=function(t,e,n){var i=this._nodes.addNode(e),r=i.getLabel();null===r?i._label=new Pe(t,n):r.setLocation(t,n)},e.prototype.createEdgeSetIntersector=function(){return new pr},e.prototype.addSelfIntersectionNodes=function(t){for(var e=this._edges.iterator();e.hasNext();)for(var n=e.next(),i=n.getLabel().getLocation(t),r=n.eiList.iterator();r.hasNext();){var o=r.next();this.addSelfIntersectionNode(t,o.coord,i)}},e.prototype.add=function(){if(1!==arguments.length)return t.prototype.add.apply(this,arguments);var e=arguments[0];if(e.isEmpty())return null;if(e instanceof ne&&(this._useBoundaryDeterminationRule=!1),e instanceof $t)this.addPolygon(e);else if(e instanceof Kt)this.addLineString(e);else if(e instanceof Qt)this.addPoint(e);else if(e instanceof te)this.addCollection(e);else if(e instanceof Xt)this.addCollection(e);else if(e instanceof ne)this.addCollection(e);else{if(!(e instanceof zt))throw new Error(e.getClass().getName());this.addCollection(e)}},e.prototype.addCollection=function(t){for(var e=0;e<t.getNumGeometries();e++){var n=t.getGeometryN(e);this.add(n)}},e.prototype.locate=function(t){return T(this._parentGeom,Zt)&&this._parentGeom.getNumGeometries()>50?(null===this._areaPtLocator&&(this._areaPtLocator=new vr(this._parentGeom)),this._areaPtLocator.locate(t)):this._ptLocator.locate(t,this._parentGeom)},e.prototype.findEdge=function(){if(1===arguments.length){var e=arguments[0];return this._lineEdgeMap.get(e)}return t.prototype.findEdge.apply(this,arguments)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e.determineBoundary=function(t,e){return t.isInBoundary(e)?w.BOUNDARY:w.INTERIOR},e}(Ye),Cr=function(){if(this._li=new rt,this._resultPrecisionModel=null,this._arg=null,1===arguments.length){var t=arguments[0];this.setComputationPrecision(t.getPrecisionModel()),this._arg=new Array(1).fill(null),this._arg[0]=new Nr(0,t)}else if(2===arguments.length){var e=arguments[0],n=arguments[1],i=gt.OGC_SFS_BOUNDARY_RULE;e.getPrecisionModel().compareTo(n.getPrecisionModel())>=0?this.setComputationPrecision(e.getPrecisionModel()):this.setComputationPrecision(n.getPrecisionModel()),this._arg=new Array(2).fill(null),this._arg[0]=new Nr(0,e,i),this._arg[1]=new Nr(1,n,i)}else if(3===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2];r.getPrecisionModel().compareTo(o.getPrecisionModel())>=0?this.setComputationPrecision(r.getPrecisionModel()):this.setComputationPrecision(o.getPrecisionModel()),this._arg=new Array(2).fill(null),this._arg[0]=new Nr(0,r,s),this._arg[1]=new Nr(1,o,s)}};Cr.prototype.getArgGeometry=function(t){return this._arg[t].getGeometry()},Cr.prototype.setComputationPrecision=function(t){this._resultPrecisionModel=t,this._li.setPrecisionModel(this._resultPrecisionModel)},Cr.prototype.interfaces_=function(){return[]},Cr.prototype.getClass=function(){return Cr};var Sr=function(){};Sr.prototype.interfaces_=function(){return[]},Sr.prototype.getClass=function(){return Sr},Sr.map=function(){if(arguments[0]instanceof ct&&T(arguments[1],Sr.MapOp)){for(var t=arguments[0],e=arguments[1],n=new Nt,i=0;i<t.getNumGeometries();i++){var r=e.map(t.getGeometryN(i));null!==r&&n.add(r)}return t.getFactory().buildGeometry(n)}if(T(arguments[0],It)&&T(arguments[1],Sr.MapOp)){for(var o=arguments[0],s=arguments[1],a=new Nt,u=o.iterator();u.hasNext();){var l=u.next(),c=s.map(l);null!==c&&a.add(c)}return a}},Sr.MapOp=function(){};var Lr=function(t){function e(){var e=arguments[0],n=arguments[1];t.call(this,e,n),this._ptLocator=new Si,this._geomFact=null,this._resultGeom=null,this._graph=null,this._edgeList=new Hn,this._resultPolyList=new Nt,this._resultLineList=new Nt,this._resultPointList=new Nt,this._graph=new Ye(new kn),this._geomFact=e.getFactory()}return t&&(e.__proto__=t),e.prototype=Object.create(t&&t.prototype),e.prototype.constructor=e,e.prototype.insertUniqueEdge=function(t){var e=this._edgeList.findEqualEdge(t);if(null!==e){var n=e.getLabel(),i=t.getLabel();e.isPointwiseEqual(t)||(i=new Pe(t.getLabel())).flip();var r=e.getDepth();r.isNull()&&r.add(n),r.add(i),n.merge(i)}else this._edgeList.add(t)},e.prototype.getGraph=function(){return this._graph},e.prototype.cancelDuplicateResultEdges=function(){for(var t=this._graph.getEdgeEnds().iterator();t.hasNext();){var e=t.next(),n=e.getSym();e.isInResult()&&n.isInResult()&&(e.setInResult(!1),n.setInResult(!1))}},e.prototype.isCoveredByLA=function(t){return!!this.isCovered(t,this._resultLineList)||!!this.isCovered(t,this._resultPolyList)},e.prototype.computeGeometry=function(t,n,i,r){var o=new Nt;return o.addAll(t),o.addAll(n),o.addAll(i),o.isEmpty()?e.createEmptyResult(r,this._arg[0].getGeometry(),this._arg[1].getGeometry(),this._geomFact):this._geomFact.buildGeometry(o)},e.prototype.mergeSymLabels=function(){for(var t=this._graph.getNodes().iterator();t.hasNext();){t.next().getEdges().mergeSymLabels()}},e.prototype.isCovered=function(t,e){for(var n=e.iterator();n.hasNext();){var i=n.next();if(this._ptLocator.locate(t,i)!==w.EXTERIOR)return!0}return!1},e.prototype.replaceCollapsedEdges=function(){for(var t=new Nt,e=this._edgeList.iterator();e.hasNext();){var n=e.next();n.isCollapsed()&&(e.remove(),t.add(n.getCollapsedEdge()))}this._edgeList.addAll(t)},e.prototype.updateNodeLabelling=function(){for(var t=this._graph.getNodes().iterator();t.hasNext();){var e=t.next(),n=e.getEdges().getLabel();e.getLabel().merge(n)}},e.prototype.getResultGeometry=function(t){return this.computeOverlay(t),this._resultGeom},e.prototype.insertUniqueEdges=function(t){for(var e=t.iterator();e.hasNext();){var n=e.next();this.insertUniqueEdge(n)}},e.prototype.computeOverlay=function(t){this.copyPoints(0),this.copyPoints(1),this._arg[0].computeSelfNodes(this._li,!1),this._arg[1].computeSelfNodes(this._li,!1),this._arg[0].computeEdgeIntersections(this._arg[1],this._li,!0);var e=new Nt;this._arg[0].computeSplitEdges(e),this._arg[1].computeSplitEdges(e),this.insertUniqueEdges(e),this.computeLabelsFromDepths(),this.replaceCollapsedEdges(),Yi.checkValid(this._edgeList.getEdges()),this._graph.addEdges(this._edgeList.getEdges()),this.computeLabelling(),this.labelIncompleteNodes(),this.findResultAreaEdges(t),this.cancelDuplicateResultEdges();var n=new ke(this._geomFact);n.add(this._graph),this._resultPolyList=n.getPolygons();var i=new ji(this,this._geomFact,this._ptLocator);this._resultLineList=i.build(t);var r=new Hi(this,this._geomFact,this._ptLocator);this._resultPointList=r.build(t),this._resultGeom=this.computeGeometry(this._resultPointList,this._resultLineList,this._resultPolyList,t)},e.prototype.labelIncompleteNode=function(t,e){var n=this._ptLocator.locate(t.getCoordinate(),this._arg[e].getGeometry());t.getLabel().setLocation(e,n)},e.prototype.copyPoints=function(t){for(var e=this._arg[t].getNodeIterator();e.hasNext();){var n=e.next();this._graph.addNode(n.getCoordinate()).setLabel(t,n.getLabel().getLocation(t))}},e.prototype.findResultAreaEdges=function(t){for(var n=this._graph.getEdgeEnds().iterator();n.hasNext();){var i=n.next(),r=i.getLabel();r.isArea()&&!i.isInteriorAreaEdge()&&e.isResultOfOp(r.getLocation(0,Se.RIGHT),r.getLocation(1,Se.RIGHT),t)&&i.setInResult(!0)}},e.prototype.computeLabelsFromDepths=function(){for(var t=this._edgeList.iterator();t.hasNext();){var e=t.next(),n=e.getLabel(),i=e.getDepth();if(!i.isNull()){i.normalize();for(var r=0;r<2;r++)n.isNull(r)||!n.isArea()||i.isNull(r)||(0===i.getDelta(r)?n.toLine(r):(et.isTrue(!i.isNull(r,Se.LEFT),"depth of LEFT side has not been initialized"),n.setLocation(r,Se.LEFT,i.getLocation(r,Se.LEFT)),et.isTrue(!i.isNull(r,Se.RIGHT),"depth of RIGHT side has not been initialized"),n.setLocation(r,Se.RIGHT,i.getLocation(r,Se.RIGHT))))}}},e.prototype.computeLabelling=function(){for(var t=this._graph.getNodes().iterator();t.hasNext();){t.next().getEdges().computeLabelling(this._arg)}this.mergeSymLabels(),this.updateNodeLabelling()},e.prototype.labelIncompleteNodes=function(){for(var t=this._graph.getNodes().iterator();t.hasNext();){var e=t.next(),n=e.getLabel();e.isIsolated()&&(n.isNull(0)?this.labelIncompleteNode(e,0):this.labelIncompleteNode(e,1)),e.getEdges().updateLabelling(n)}},e.prototype.isCoveredByA=function(t){return!!this.isCovered(t,this._resultPolyList)},e.prototype.interfaces_=function(){return[]},e.prototype.getClass=function(){return e},e}(Cr);Lr.overlayOp=function(t,e,n){return new Lr(t,e).getResultGeometry(n)},Lr.intersection=function(t,e){if(t.isEmpty()||e.isEmpty())return Lr.createEmptyResult(Lr.INTERSECTION,t,e,t.getFactory());if(t.isGeometryCollection()){var n=e;return ki.map(t,{interfaces_:function(){return[Sr.MapOp]},map:function(t){return t.intersection(n)}})}return t.checkNotGeometryCollection(t),t.checkNotGeometryCollection(e),or.overlayOp(t,e,Lr.INTERSECTION)},Lr.symDifference=function(t,e){if(t.isEmpty()||e.isEmpty()){if(t.isEmpty()&&e.isEmpty())return Lr.createEmptyResult(Lr.SYMDIFFERENCE,t,e,t.getFactory());if(t.isEmpty())return e.copy();if(e.isEmpty())return t.copy()}return t.checkNotGeometryCollection(t),t.checkNotGeometryCollection(e),or.overlayOp(t,e,Lr.SYMDIFFERENCE)},Lr.resultDimension=function(t,e,n){var i=e.getDimension(),r=n.getDimension(),o=-1;switch(t){case Lr.INTERSECTION:o=Math.min(i,r);break;case Lr.UNION:o=Math.max(i,r);break;case Lr.DIFFERENCE:o=i;break;case Lr.SYMDIFFERENCE:o=Math.max(i,r)}return o},Lr.createEmptyResult=function(t,e,n,i){var r=null;switch(Lr.resultDimension(t,e,n)){case-1:r=i.createGeometryCollection(new Array(0).fill(null));break;case 0:r=i.createPoint();break;case 1:r=i.createLineString();break;case 2:r=i.createPolygon()}return r},Lr.difference=function(t,e){return t.isEmpty()?Lr.createEmptyResult(Lr.DIFFERENCE,t,e,t.getFactory()):e.isEmpty()?t.copy():(t.checkNotGeometryCollection(t),t.checkNotGeometryCollection(e),or.overlayOp(t,e,Lr.DIFFERENCE))},Lr.isResultOfOp=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1],n=t.getLocation(0),i=t.getLocation(1);return Lr.isResultOfOp(n,i,e)}if(3===arguments.length){var r=arguments[0],o=arguments[1],s=arguments[2];switch(r===w.BOUNDARY&&(r=w.INTERIOR),o===w.BOUNDARY&&(o=w.INTERIOR),s){case Lr.INTERSECTION:return r===w.INTERIOR&&o===w.INTERIOR;case Lr.UNION:return r===w.INTERIOR||o===w.INTERIOR;case Lr.DIFFERENCE:return r===w.INTERIOR&&o!==w.INTERIOR;case Lr.SYMDIFFERENCE:return r===w.INTERIOR&&o!==w.INTERIOR||r!==w.INTERIOR&&o===w.INTERIOR}return!1}},Lr.INTERSECTION=1,Lr.UNION=2,Lr.DIFFERENCE=3,Lr.SYMDIFFERENCE=4;var br=function(){this._g=null,this._boundaryDistanceTolerance=null,this._linework=null,this._ptLocator=new Si,this._seg=new dn;var t=arguments[0],e=arguments[1];this._g=t,this._boundaryDistanceTolerance=e,this._linework=this.extractLinework(t)};br.prototype.isWithinToleranceOfBoundary=function(t){for(var e=0;e<this._linework.getNumGeometries();e++)for(var n=this._linework.getGeometryN(e).getCoordinateSequence(),i=0;i<n.size()-1;i++){n.getCoordinate(i,this._seg.p0),n.getCoordinate(i+1,this._seg.p1);if(this._seg.distance(t)<=this._boundaryDistanceTolerance)return!0}return!1},br.prototype.getLocation=function(t){return this.isWithinToleranceOfBoundary(t)?w.BOUNDARY:this._ptLocator.locate(t,this._g)},br.prototype.extractLinework=function(t){var e=new wr;t.apply(e);var n=e.getLinework(),i=_e.toLineStringArray(n);return t.getFactory().createMultiLineString(i)},br.prototype.interfaces_=function(){return[]},br.prototype.getClass=function(){return br};var wr=function(){this._linework=null,this._linework=new Nt};wr.prototype.getLinework=function(){return this._linework},wr.prototype.filter=function(t){if(t instanceof $t){var e=t;this._linework.add(e.getExteriorRing());for(var n=0;n<e.getNumInteriorRing();n++)this._linework.add(e.getInteriorRingN(n))}},wr.prototype.interfaces_=function(){return[Vt]},wr.prototype.getClass=function(){return wr};var Or=function(){this._g=null,this._doLeft=!0,this._doRight=!0;var t=arguments[0];this._g=t};Or.prototype.extractPoints=function(t,e,n){for(var i=t.getCoordinates(),r=0;r<i.length-1;r++)this.computeOffsetPoints(i[r],i[r+1],e,n)},Or.prototype.setSidesToGenerate=function(t,e){this._doLeft=t,this._doRight=e},Or.prototype.getPoints=function(t){for(var e=new Nt,n=Ci.getLines(this._g).iterator();n.hasNext();){var i=n.next();this.extractPoints(i,t,e)}return e},Or.prototype.computeOffsetPoints=function(t,e,n,i){var r=e.x-t.x,o=e.y-t.y,s=Math.sqrt(r*r+o*o),a=n*r/s,u=n*o/s,l=(e.x+t.x)/2,c=(e.y+t.y)/2;if(this._doLeft){var p=new C(l-u,c+a);i.add(p)}if(this._doRight){var h=new C(l+u,c-a);i.add(h)}},Or.prototype.interfaces_=function(){return[]},Or.prototype.getClass=function(){return Or};var Tr=function t(){this._geom=null,this._locFinder=null,this._location=new Array(3).fill(null),this._invalidLocation=null,this._boundaryDistanceTolerance=t.TOLERANCE,this._testCoords=new Nt;var e=arguments[0],n=arguments[1],i=arguments[2];this._boundaryDistanceTolerance=t.computeBoundaryDistanceTolerance(e,n),this._geom=[e,n,i],this._locFinder=[new br(this._geom[0],this._boundaryDistanceTolerance),new br(this._geom[1],this._boundaryDistanceTolerance),new br(this._geom[2],this._boundaryDistanceTolerance)]},Rr={TOLERANCE:{configurable:!0}};Tr.prototype.reportResult=function(t,e,n){Y.out.println("Overlay result invalid - A:"+w.toLocationSymbol(e[0])+" B:"+w.toLocationSymbol(e[1])+" expected:"+(n?"i":"e")+" actual:"+w.toLocationSymbol(e[2]))},Tr.prototype.isValid=function(t){this.addTestPts(this._geom[0]),this.addTestPts(this._geom[1]);var e=this.checkValid(t);return e},Tr.prototype.checkValid=function(){if(1===arguments.length){for(var t=arguments[0],e=0;e<this._testCoords.size();e++){var n=this._testCoords.get(e);if(!this.checkValid(t,n))return this._invalidLocation=n,!1}return!0}if(2===arguments.length){var i=arguments[0],r=arguments[1];return this._location[0]=this._locFinder[0].getLocation(r),this._location[1]=this._locFinder[1].getLocation(r),this._location[2]=this._locFinder[2].getLocation(r),!!Tr.hasLocation(this._location,w.BOUNDARY)||this.isValidResult(i,this._location)}},Tr.prototype.addTestPts=function(t){var e=new Or(t);this._testCoords.addAll(e.getPoints(5*this._boundaryDistanceTolerance))},Tr.prototype.isValidResult=function(t,e){var n=Lr.isResultOfOp(e[0],e[1],t),i=!(n^e[2]===w.INTERIOR);return i||this.reportResult(t,e,n),i},Tr.prototype.getInvalidLocation=function(){return this._invalidLocation},Tr.prototype.interfaces_=function(){return[]},Tr.prototype.getClass=function(){return Tr},Tr.hasLocation=function(t,e){for(var n=0;n<3;n++)if(t[n]===e)return!0;return!1},Tr.computeBoundaryDistanceTolerance=function(t,e){return Math.min(Ji.computeSizeBasedSnapTolerance(t),Ji.computeSizeBasedSnapTolerance(e))},Tr.isValid=function(t,e,n,i){return new Tr(t,e,i).isValid(n)},Rr.TOLERANCE.get=function(){return 1e-6},Object.defineProperties(Tr,Rr);var Pr=function t(e){this._geomFactory=null,this._skipEmpty=!1,this._inputGeoms=null,this._geomFactory=t.extractFactory(e),this._inputGeoms=e};Pr.prototype.extractElements=function(t,e){if(null===t)return null;for(var n=0;n<t.getNumGeometries();n++){var i=t.getGeometryN(n);this._skipEmpty&&i.isEmpty()||e.add(i)}},Pr.prototype.combine=function(){for(var t=new Nt,e=this._inputGeoms.iterator();e.hasNext();){var n=e.next();this.extractElements(n,t)}return 0===t.size()?null!==this._geomFactory?this._geomFactory.createGeometryCollection(null):null:this._geomFactory.buildGeometry(t)},Pr.prototype.interfaces_=function(){return[]},Pr.prototype.getClass=function(){return Pr},Pr.combine=function(){if(1===arguments.length){var t=arguments[0];return new Pr(t).combine()}if(2===arguments.length){var e=arguments[0],n=arguments[1];return new Pr(Pr.createList(e,n)).combine()}if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2];return new Pr(Pr.createList(i,r,o)).combine()}},Pr.extractFactory=function(t){return t.isEmpty()?null:t.iterator().next().getFactory()},Pr.createList=function(){if(2===arguments.length){var t=arguments[0],e=arguments[1],n=new Nt;return n.add(t),n.add(e),n}if(3===arguments.length){var i=arguments[0],r=arguments[1],o=arguments[2],s=new Nt;return s.add(i),s.add(r),s.add(o),s}};var Dr=function(){this._inputPolys=null,this._geomFactory=null;var t=arguments[0];this._inputPolys=t,null===this._inputPolys&&(this._inputPolys=new Nt)},Mr={STRTREE_NODE_CAPACITY:{configurable:!0}};Dr.prototype.reduceToGeometries=function(t){for(var e=new Nt,n=t.iterator();n.hasNext();){var i=n.next(),r=null;T(i,xt)?r=this.unionTree(i):i instanceof ct&&(r=i),e.add(r)}return e},Dr.prototype.extractByEnvelope=function(t,e,n){for(var i=new Nt,r=0;r<e.getNumGeometries();r++){var o=e.getGeometryN(r);o.getEnvelopeInternal().intersects(t)?i.add(o):n.add(o)}return this._geomFactory.buildGeometry(i)},Dr.prototype.unionOptimized=function(t,e){var n=t.getEnvelopeInternal(),i=e.getEnvelopeInternal();if(!n.intersects(i)){return Pr.combine(t,e)}if(t.getNumGeometries()<=1&&e.getNumGeometries()<=1)return this.unionActual(t,e);var r=n.intersection(i);return this.unionUsingEnvelopeIntersection(t,e,r)},Dr.prototype.union=function(){if(null===this._inputPolys)throw new Error("union() method cannot be called twice");if(this._inputPolys.isEmpty())return null;this._geomFactory=this._inputPolys.iterator().next().getFactory();for(var t=new sn(Dr.STRTREE_NODE_CAPACITY),e=this._inputPolys.iterator();e.hasNext();){var n=e.next();t.insert(n.getEnvelopeInternal(),n)}this._inputPolys=null;var i=t.itemsTree();return this.unionTree(i)},Dr.prototype.binaryUnion=function(){if(1===arguments.length){var t=arguments[0];return this.binaryUnion(t,0,t.size())}if(3===arguments.length){var e=arguments[0],n=arguments[1],i=arguments[2];if(i-n<=1){var r=Dr.getGeometry(e,n);return this.unionSafe(r,null)}if(i-n==2)return this.unionSafe(Dr.getGeometry(e,n),Dr.getGeometry(e,n+1));var o=Math.trunc((i+n)/2),s=this.binaryUnion(e,n,o),a=this.binaryUnion(e,o,i);return this.unionSafe(s,a)}},Dr.prototype.repeatedUnion=function(t){for(var e=null,n=t.iterator();n.hasNext();){var i=n.next();e=null===e?i.copy():e.union(i)}return e},Dr.prototype.unionSafe=function(t,e){return null===t&&null===e?null:null===t?e.copy():null===e?t.copy():this.unionOptimized(t,e)},Dr.prototype.unionActual=function(t,e){return Dr.restrictToPolygons(t.union(e))},Dr.prototype.unionTree=function(t){var e=this.reduceToGeometries(t);return this.binaryUnion(e)},Dr.prototype.unionUsingEnvelopeIntersection=function(t,e,n){var i=new Nt,r=this.extractByEnvelope(n,t,i),o=this.extractByEnvelope(n,e,i),s=this.unionActual(r,o);i.add(s);return Pr.combine(i)},Dr.prototype.bufferUnion=function(){if(1===arguments.length){var t=arguments[0];return t.get(0).getFactory().buildGeometry(t).buffer(0)}if(2===arguments.length){var e=arguments[0],n=arguments[1];return e.getFactory().createGeometryCollection([e,n]).buffer(0)}},Dr.prototype.interfaces_=function(){return[]},Dr.prototype.getClass=function(){return Dr},Dr.restrictToPolygons=function(t){if(T(t,Zt))return t;var e=Ni.getPolygons(t);return 1===e.size()?e.get(0):t.getFactory().createMultiPolygon(_e.toPolygonArray(e))},Dr.getGeometry=function(t,e){return e>=t.size()?null:t.get(e)},Dr.union=function(t){return new Dr(t).union()},Mr.STRTREE_NODE_CAPACITY.get=function(){return 4},Object.defineProperties(Dr,Mr);var Ar=function(){};Ar.prototype.interfaces_=function(){return[]},Ar.prototype.getClass=function(){return Ar},Ar.union=function(t,e){if(t.isEmpty()||e.isEmpty()){if(t.isEmpty()&&e.isEmpty())return Lr.createEmptyResult(Lr.UNION,t,e,t.getFactory());if(t.isEmpty())return e.copy();if(e.isEmpty())return t.copy()}return t.checkNotGeometryCollection(t),t.checkNotGeometryCollection(e),or.overlayOp(t,e,Lr.UNION)},t.GeoJSONReader=Ne,t.GeoJSONWriter=Ce,t.OverlayOp=Lr,t.UnionOp=Ar,t.BufferOp=di,Object.defineProperty(t,"__esModule",{value:!0})});

},{}]},{},[1]);
