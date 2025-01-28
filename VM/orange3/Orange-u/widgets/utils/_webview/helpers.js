/**
 * Our general helpers for Highcharts, JS, QWebView bridge ...
 */

// Prevent descent into the values of these keys
var _PREVENT_KEYS = [
    'data'  // This is the numeric payload in Highcharts
];

function fixupPythonObject(obj) {
    /**
     * Replace any strings with their eval'd value if they
     * start with an empty JavaScript block comment, i.e. these
     * four characters:
     */                 /**/
    if (typeof obj === 'undefined' || obj === null)
        return;

    var keys = Object.keys(obj);
    for (var i=0; i<keys.length; ++i) {
        var key = keys[i],
            val = obj[key];

        if (val === null || val === undefined)
            continue;

        // Make sure arrays are of type Array and not Qt's RuntimeArray
        // Can probably be removed once Qt 4's WebKit support is dropped.
        if (val.constructor === Array &&
            !(Object.prototype.toString.call(val) == '[object Array]'))
            // FIXME: This is suboptimal, but what can we do? Simple "casting"
            // into Array with Array.prototype.slice() didn't seem to work,
            // and using Array.prototype.map is much slower.
            obj[key] = val = JSON.parse(JSON.stringify(val));

        if (typeof val === 'string' && val.indexOf('/**/') == 0) {
            obj[key] = eval(val)
        } else if (val.constructor === Object ||
                   val.constructor === Array && _PREVENT_KEYS.indexOf(key) == -1)
            fixupPythonObject(val);
    }
}


window.addEventListener("keyup", function(ev) {
    var ESCAPE_KEY = 27;
    // Cmd+W on OS X seems to work OK without this
    if (ev.keyCode == ESCAPE_KEY ||
        ev.code == 'Escape' || ev.key == 'Escape')
        window.__bridge.hideWindow();
});
