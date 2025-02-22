//  ordered_queue_map.h
//
//  Copyright 2024 Antoine Passemiers <antoine.passemiers@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA 02110-1301, USA.

#ifndef MODALIB_ORDERED_QUEUE_MAP_H
#define MODALIB_ORDERED_QUEUE_MAP_H

#include <cstdint>
#include <queue>
#include <unordered_map>
#include <iostream>

namespace modalib {

template<class K, class V>
class OrderedQueueMap {
private:
    size_t max_size;
    std::queue<K> keys;
    std::unordered_map<K, V> values;
public:
    OrderedQueueMap() = delete;

    explicit OrderedQueueMap(size_t max_size) : max_size(max_size), keys(), values() {
        values.reserve(max_size);
    }

    ~OrderedQueueMap() = default;

    void add(const K &key, const V &value) {
        if (!contains(key)) {
            this->keys.push(key);
            this->values[key] = value;
        }
        this->update();
    }

    std::pair<K, V> front() const {
        const K &key = this->keys.front();
        const V &value = this->values.at(key);
        return {key, value};
    }

private:
    void update() {
        if (this->keys.size() > 2 * this->max_size) {
            std::queue<K> new_queue;
            while (!this->keys.empty()) {
                const K& key = this->keys.front();
                if (this->contains_value(key)) {
                    new_queue.push(key);
                }
                this->keys.pop();
            }
            this->keys = new_queue;
        }
        while ((this->keys.size() > 0) && (!this->contains_value(this->keys.front()))) {
            this->keys.pop();
        }
    }

    inline bool contains_value(const K &key) const { return this->values.find(key) != this->values.end(); }

public:
    inline void remove(const K &key) {
        this->values.erase(key);
        update();
    }

    inline bool contains(const K &key) const { return this->contains_value(key); }

    inline const V& get(const K &key) const { return this->values.at(key); }

    inline size_t size() const { return this->values.size(); }

    inline bool size_exceeded() const { return this->values.size() >= this->max_size; }

    void pop() {
        if (this->keys.size() == 0) {
            this->values.clear();
        } else {
            const K &key = this->keys.front();
            this->values.erase(key);
            this->keys.pop();
        }
    }

    void print_debug_info() const {
        std::cout << keys.size() << ", " << values.size() << std::endl;
    }
};

}  // namespace modalib

#endif //MODALIB_ORDERED_QUEUE_MAP_H
