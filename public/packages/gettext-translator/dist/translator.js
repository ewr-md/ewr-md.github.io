(function (global, factory) {
  if (typeof define === "function" && define.amd) {
    define(["exports"], factory);
  } else if (typeof exports !== "undefined") {
    factory(exports);
  } else {
    var mod = {
      exports: {}
    };
    factory(mod.exports);
    global.translator = mod.exports;
  }
})(typeof globalThis !== "undefined" ? globalThis : typeof self !== "undefined" ? self : this, function (_exports) {
  "use strict";

  Object.defineProperty(_exports, "__esModule", {
    value: true
  });
  _exports["default"] = void 0;

  function _typeof(obj) { if (typeof Symbol === "function" && typeof Symbol.iterator === "symbol") { _typeof = function _typeof(obj) { return typeof obj; }; } else { _typeof = function _typeof(obj) { return obj && typeof Symbol === "function" && obj.constructor === Symbol && obj !== Symbol.prototype ? "symbol" : typeof obj; }; } return _typeof(obj); }

  function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

  function _defineProperties(target, props) { for (var i = 0; i < props.length; i++) { var descriptor = props[i]; descriptor.enumerable = descriptor.enumerable || false; descriptor.configurable = true; if ("value" in descriptor) descriptor.writable = true; Object.defineProperty(target, descriptor.key, descriptor); } }

  function _createClass(Constructor, protoProps, staticProps) { if (protoProps) _defineProperties(Constructor.prototype, protoProps); if (staticProps) _defineProperties(Constructor, staticProps); return Constructor; }

  var Translator =
  /*#__PURE__*/
  function () {
    function Translator(translations) {
      _classCallCheck(this, Translator);

      this.dictionary = {};
      this.plurals = {};
      this.domain = undefined;

      if (translations) {
        this.loadTranslations(translations);
      }
    }

    _createClass(Translator, [{
      key: "loadTranslations",
      value: function loadTranslations(translations) {
        var domain = translations.domain || '';

        if (this.domain === undefined) {
          this.domain = domain;
        }

        if (this.dictionary[domain]) {
          mergeTranslations(this.dictionary[domain], translations.messages);
          return this;
        }

        if (translations.fn) {
          this.plurals[domain] = {
            fn: translations.fn
          };
        } else if (translations['plural-forms']) {
          var plural = translations['plural-forms'].split(';', 2);
          this.plurals[domain] = {
            count: parseInt(plural[0].replace('nplurals=', '')),
            code: plural[1].replace('plural=', 'return ') + ';'
          };
        }

        this.dictionary[domain] = translations.messages;
        return this;
      }
    }, {
      key: "defaultDomain",
      value: function defaultDomain(domain) {
        this.domain = domain;
        return this;
      }
    }, {
      key: "gettext",
      value: function gettext(original) {
        for (var _len = arguments.length, args = new Array(_len > 1 ? _len - 1 : 0), _key = 1; _key < _len; _key++) {
          args[_key - 1] = arguments[_key];
        }

        return this.format.apply(this, [this.translate(undefined, undefined, original)].concat(args));
      }
    }, {
      key: "ngettext",
      value: function ngettext(original, plural, counter) {
        for (var _len2 = arguments.length, args = new Array(_len2 > 3 ? _len2 - 3 : 0), _key2 = 3; _key2 < _len2; _key2++) {
          args[_key2 - 3] = arguments[_key2];
        }

        return this.format.apply(this, [this.translatePlural(undefined, undefined, original, plural, counter)].concat(args));
      }
    }, {
      key: "dngettext",
      value: function dngettext(domain, original, plural, counter) {
        for (var _len3 = arguments.length, args = new Array(_len3 > 4 ? _len3 - 4 : 0), _key3 = 4; _key3 < _len3; _key3++) {
          args[_key3 - 4] = arguments[_key3];
        }

        return this.format.apply(this, [this.translatePlural(domain, undefined, original, plural, counter)].concat(args));
      }
    }, {
      key: "npgettext",
      value: function npgettext(context, original, plural, counter) {
        for (var _len4 = arguments.length, args = new Array(_len4 > 4 ? _len4 - 4 : 0), _key4 = 4; _key4 < _len4; _key4++) {
          args[_key4 - 4] = arguments[_key4];
        }

        return this.format.apply(this, [this.translatePlural(undefined, context, original, plural, counter)].concat(args));
      }
    }, {
      key: "pgettext",
      value: function pgettext(context, original) {
        for (var _len5 = arguments.length, args = new Array(_len5 > 2 ? _len5 - 2 : 0), _key5 = 2; _key5 < _len5; _key5++) {
          args[_key5 - 2] = arguments[_key5];
        }

        return this.format.apply(this, [this.translate(undefined, context, original)].concat(args));
      }
    }, {
      key: "dgettext",
      value: function dgettext(domain, original) {
        for (var _len6 = arguments.length, args = new Array(_len6 > 2 ? _len6 - 2 : 0), _key6 = 2; _key6 < _len6; _key6++) {
          args[_key6 - 2] = arguments[_key6];
        }

        return this.format.apply(this, [this.translate(domain, undefined, original)].concat(args));
      }
    }, {
      key: "dpgettext",
      value: function dpgettext(domain, context, original) {
        for (var _len7 = arguments.length, args = new Array(_len7 > 3 ? _len7 - 3 : 0), _key7 = 3; _key7 < _len7; _key7++) {
          args[_key7 - 3] = arguments[_key7];
        }

        return this.format.apply(this, [this.translate(domain, context, original)].concat(args));
      }
    }, {
      key: "dnpgettext",
      value: function dnpgettext(domain, context, original, plural, counter) {
        for (var _len8 = arguments.length, args = new Array(_len8 > 5 ? _len8 - 5 : 0), _key8 = 5; _key8 < _len8; _key8++) {
          args[_key8 - 5] = arguments[_key8];
        }

        return this.format.apply(this, [this.translatePlural(domain, context, original, plural, counter)].concat(args));
      }
    }, {
      key: "__",
      value: function __(original) {
        for (var _len9 = arguments.length, args = new Array(_len9 > 1 ? _len9 - 1 : 0), _key9 = 1; _key9 < _len9; _key9++) {
          args[_key9 - 1] = arguments[_key9];
        }

        return this.gettext.apply(this, [original].concat(args));
      }
    }, {
      key: "n__",
      value: function n__(original, plural, value) {
        for (var _len10 = arguments.length, args = new Array(_len10 > 3 ? _len10 - 3 : 0), _key10 = 3; _key10 < _len10; _key10++) {
          args[_key10 - 3] = arguments[_key10];
        }

        return this.ngettext.apply(this, [original, plural, value].concat(args));
      }
    }, {
      key: "p__",
      value: function p__(context, original) {
        for (var _len11 = arguments.length, args = new Array(_len11 > 2 ? _len11 - 2 : 0), _key11 = 2; _key11 < _len11; _key11++) {
          args[_key11 - 2] = arguments[_key11];
        }

        return this.pgettext.apply(this, [context, original].concat(args));
      }
    }, {
      key: "d__",
      value: function d__(domain, original) {
        for (var _len12 = arguments.length, args = new Array(_len12 > 2 ? _len12 - 2 : 0), _key12 = 2; _key12 < _len12; _key12++) {
          args[_key12 - 2] = arguments[_key12];
        }

        return this.dgettext.apply(this, [domain, original].concat(args));
      }
    }, {
      key: "dp__",
      value: function dp__(domain, context, original) {
        for (var _len13 = arguments.length, args = new Array(_len13 > 3 ? _len13 - 3 : 0), _key13 = 3; _key13 < _len13; _key13++) {
          args[_key13 - 3] = arguments[_key13];
        }

        return this.dpgettext.apply(this, [domain, context, original].concat(args));
      }
    }, {
      key: "np__",
      value: function np__(context, original, plural, value) {
        for (var _len14 = arguments.length, args = new Array(_len14 > 4 ? _len14 - 4 : 0), _key14 = 4; _key14 < _len14; _key14++) {
          args[_key14 - 4] = arguments[_key14];
        }

        return this.npgettext.apply(this, [context, original, plural, value].concat(args));
      }
    }, {
      key: "dnp__",
      value: function dnp__(domain, context, original, plural, value) {
        for (var _len15 = arguments.length, args = new Array(_len15 > 5 ? _len15 - 5 : 0), _key15 = 5; _key15 < _len15; _key15++) {
          args[_key15 - 5] = arguments[_key15];
        }

        return this.dnpgettext.apply(this, [domain, context, original, plural, value].concat(args));
      }
    }, {
      key: "format",
      value: function format(text) {
        for (var _len16 = arguments.length, args = new Array(_len16 > 1 ? _len16 - 1 : 0), _key16 = 1; _key16 < _len16; _key16++) {
          args[_key16 - 1] = arguments[_key16];
        }

        if (!args.length) {
          return text;
        }

        if (_typeof(args[0]) === 'object') {
          Object.keys(args[0]).forEach(function (search) {
            text = text.replace(search, args[0][search]);
          });
          return text;
        }

        return text.replace(/(%[sd])/g, function (match) {
          if (!args.length) {
            return match;
          }

          switch (match) {
            case '%s':
              return args.shift();

            case '%d':
              return parseFloat(args.shift());
          }
        });
      }
    }, {
      key: "translate",
      value: function translate(domain, context, original) {
        var translation = this.getTranslation(domain, context, original);
        return translation && translation[0] ? translation[0] : original;
      }
    }, {
      key: "translatePlural",
      value: function translatePlural(domain, context, original, plural, counter) {
        var translation = this.getTranslation(domain, context, original);
        var index = this.getPluralIndex(domain, counter);
        return translation && translation[index] ? translation[index] : index === 0 ? original : plural;
      }
    }, {
      key: "getTranslation",
      value: function getTranslation(domain, context, original) {
        domain = domain || this.domain;
        context = context || '';

        if (!this.dictionary[domain] || !this.dictionary[domain][context] || !this.dictionary[domain][context][original]) {
          return undefined;
        }

        var translation = this.dictionary[domain][context][original];
        return Array.isArray(translation) ? translation : [translation];
      }
    }, {
      key: "getPluralIndex",
      value: function getPluralIndex(domain, value) {
        domain = domain || this.domain;

        if (!this.plurals[domain]) {
          return value == 1 ? 0 : 1;
        }

        if (!this.plurals[domain].fn) {
          this.plurals[domain].fn = new Function('n', this.plurals[domain].code);
        }

        return this.plurals[domain].fn.call(this, value) + 0;
      }
    }]);

    return Translator;
  }();

  _exports["default"] = Translator;

  function mergeTranslations(translations, newTranslations) {
    for (var context in newTranslations) {
      if (!translations[context]) {
        translations[context] = newTranslations[context];
        continue;
      }

      for (var original in newTranslations[context]) {
        translations[context][original] = newTranslations[context][original];
      }
    }
  }
});